// In-place DRAGEN demux/alignment and QC. Existing outputs control skipping; task caching is disabled.

// Nextflow wraps exceptions from helper functions; log the cause before throwing.
def fail(message) {
    log.error message
    error message
}

def quote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def validSample(sample) {
    if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || sample == 'pipeline_info' || sample.startsWith('multiqc')) {
        fail("Invalid sample directory name: ${sample}")
    }
}

def alignPanel(idx) {
    return ['U5N2_I10': 'HEME', 'I10': 'PCP', 'I8N2_N10': 'SKIP'][idx]
}

def sheetIndex(name) {
    return name.replace('SampleSheet_', '').replace('.csv', '')
}

def sampleSheets(runDir) {
    def sheets = files("${runDir}/*SampleSheet_*").sort { s -> s.name }
    if (!sheets) {
        fail("No samplesheets found in ${runDir}")
    }
    return sheets
}

// DRAGEN creates the index directory, so its incomplete marker lives beside it.
def demuxMarker(fastqRun, idx) {
    return "${fastqRun}/.${idx}.bclqc_demux_incomplete".toString()
}

def demuxProblems(idxs, fastqRun) {
    return idxs.collect { i ->
        def d = "${fastqRun}/${i}"
        if (file(demuxMarker(fastqRun, i)).exists()) {
            return "${i}: demux incomplete (${demuxMarker(fastqRun, i)} present: DRAGEN failed or was interrupted)".toString()
        }
        if (!file("${d}/Reports/fastq_list.csv").exists()) {
            return "${i}: not demultiplexed (${d}/Reports/fastq_list.csv missing)".toString()
        }
        return null
    }.findAll { x -> x }
}

def demuxRecovery(fastqRun) {
    return "Remove each failed index dir listed above AND its marker (rm -rf ${fastqRun}/<index> ${fastqRun}/.<index>.bclqc_demux_incomplete), " +
        "or the whole ${fastqRun}, then re-run the demux step " +
        "(e.g. --steps demux,align,qc); demux then processes only the missing index dirs. " +
        "Do NOT run --steps align,qc on an incomplete demux."
}

def checkDemuxComplete(runDir, fastqRun) {
    def bad = demuxProblems(sampleSheets(runDir).collect { s -> sheetIndex(s.name) }, fastqRun)
    if (bad) {
        fail("Demultiplexing is incomplete in ${fastqRun}; refusing to align:\n  ${bad.join('\n  ')}\n${demuxRecovery(fastqRun)}")
    }
}

def fastqDirs(fastqRun) {
    return files("${fastqRun}/*", type: 'dir').collect { d -> [d.name, d.toString()] }
}

def canonPath(p) {
    def f = file(p.toString()).toAbsolutePath().normalize()
    if (f.exists()) {
        return f.toRealPath().toString()
    }
    if (f.parent == null) {
        return f.toString()
    }
    def b = canonPath(f.parent.toString())
    return (b == '/' ? '' : b) + '/' + f.name
}

def forbiddenPrefix(path, prefixes) {
    def ps = [file(path).toAbsolutePath().normalize().toString(), canonPath(path)].unique()
    return prefixes.find { x ->
        def xs = [file(x).toAbsolutePath().normalize().toString(), canonPath(x)].unique()
        ps.any { q -> xs.any { y -> y == '/' || q == y || q.startsWith(y + '/') } }
    }
}

def alignJobs(idx, dir, bamsRun) {
    def panel = alignPanel(idx)
    if (panel == 'SKIP') {
        log.info "Skipping alignment for exome samples in ${dir}"
        return []
    }
    def lists = files("${dir}/**fastq_list.csv").findAll { f -> f.name == 'fastq_list.csv' }
    if (!lists) {
        fail("No fastq_list.csv found under ${dir}")
    }
    def jobs = []
    lists.each { fl ->
        def lines = fl.readLines().findAll { l -> l.trim() }
        if (!lines) {
            fail("Empty fastq list: ${fl}")
        }
        def col = lines[0].split(',').toList().indexOf('RGSM')
        if (col < 0) {
            fail("No RGSM column in ${fl}")
        }
        def rows = lines.drop(1).collect { l -> l.split(',', -1) }
        rows.each { f ->
            if (f.size() <= col || !f[col]) {
                fail("Missing RGSM value in ${fl}: ${f.join(',')}")
            }
        }
        rows.collect { f -> f[col] }.unique().each { s ->
            def out = "${bamsRun}/${s}".toString()
            validSample(s)
            jobs.add([s, panel, fl.toString(), out])
        }
    }
    return jobs
}

def allAlignJobs(dirs, bamsRun, skipExisting = true) {
    if (!dirs) {
        fail("No FASTQ index directories found to align")
    }
    def bad = dirs.findAll { d -> alignPanel(d[0]) == null }.collect { d -> d[0] }
    if (bad) {
        fail("Could not determine panel/bed file for index dir(s) ${bad} (known: I10, U5N2_I10, I8N2_N10). No samples were aligned.")
    }
    def seen = [:]
    def jobs = []
    dirs.sort { d -> d[0] }.each { d ->
        alignJobs(d[0], d[1], bamsRun).each { j ->
            if (seen.containsKey(j[0])) {
                fail("Sample ${j[0]} appears in both ${seen[j[0]]} and ${d[0]} fastq lists")
            }
            seen[j[0]] = d[0]
            jobs.add(j)
        }
    }
    if (!skipExisting) { return jobs }
    return jobs.findAll { j ->
        def (sample, panel, fastqs, out) = j
        if (file("${out}/.bclqc_align_incomplete").exists()) {
            fail("Partial alignment in ${out}. ${alignRecovery()}")
        }
        if (file(out).exists()) {
            if (!file("${out}/${sample}.cram").exists()) {
                fail("Alignment directory ${out} exists without ${sample}.cram; remove it before rerunning alignment")
            }
            log.info "Skipping existing alignment: ${sample}"
            return false
        }
        return true
    }
}

def alignRecovery() {
    return "Recovery: delete each marked sample dir (rm -rf <that dir>) and re-run --steps align,qc for its run. " +
        "Never run --steps qc alone after an ALIGN failure."
}

def sampleRows(rows, bamsRun) {
    if (!rows) { fail('sampleinfo has no samples') }
    def seen = [] as Set
    def out = []
    rows.eachWithIndex { row, i ->
        def sample = row['Samples']
        def panel = row['Panel'] ?: row['Tumor']
        def bam = row['BAM Path']
        if (!sample || !panel || !bam) {
            fail("sampleinfo row ${i + 2} is missing Samples / Panel (or Tumor) / BAM Path: ${row}")
        }
        validSample(sample)
        if (!seen.add(sample)) { fail("Duplicate sampleinfo sample: ${sample}") }
        out.add([i, sample, panel, file(bam.toString()).toAbsolutePath().normalize().toString(), "${bamsRun}/${sample}".toString()])
    }
    return out
}

def checkAlignComplete(rows) {
    rows.each { r ->
        if (!file(r[3]).exists()) { fail("Alignment file not found: ${r[3]}") }
    }
    def bad = rows.collect { r ->
        def ms = ["${r[4]}/.bclqc_align_incomplete", "${file(r[3]).parent}/.bclqc_align_incomplete"]*.toString().unique().findAll { m -> file(m).exists() }
        ms ? "${r[1]}: alignment incomplete (${ms.join(', ')} present: DRAGEN ALIGN failed or was interrupted)".toString() : null
    }.findAll { x -> x }
    if (bad) {
        fail("Refusing to QC partially aligned sample(s):\n  ${bad.join('\n  ')}\n${alignRecovery()}")
    }
}

// Coordinate DRAGEN across launches; read-only opening permits locks created by another user.
def dragenLock() {
    def l = params.dragen_lock
    def w = params.dragen_lock_wait
    return """
    mkdir -p "\$(dirname ${quote(l)})"
    [ -e ${quote(l)} ] || ( umask 000; : >> ${quote(l)} ) 2>/dev/null || true
    exec 8<${quote(l)} || { echo ${quote("ERROR: cannot open DRAGEN lock ${l}")} >&2; exit 1; }
    echo ${quote("Waiting for DRAGEN lock ${l} (up to ${w} s)")} >&2
    LOCK_T0=\$(date +%s)
    flock -w ${w} 8 || { echo ${quote("DRAGEN lock ${l} timed out after ${w} s; inspect with lslocks")} >&2; exit 1; }
    echo "DRAGEN lock acquired; waited \$(( \$(date +%s) - LOCK_T0 )) s" >&2
    """
}

def writeGuard(paths) {
    return params.forbid_prefixes ? "check_write_paths.py ${quote(params.forbid_prefixes)} " + paths.collect { quote(it) }.join(' ') : ''
}

def readCommand(args) {
    def process = new ProcessBuilder(args.collect { it.toString() }).redirectErrorStream(true).start()
    def output = process.inputStream.text.trim()
    if (process.waitFor() != 0) { fail(output) }
    return output
}

def dryRun(runDir, fastqRun, bamsRun) {
    def sheets = sampleSheets(runDir)
    checkDemuxComplete(runDir, fastqRun)
    def jobs = allAlignJobs(fastqDirs(fastqRun), bamsRun, false)
    def rows = file(params.sampleinfo).splitCsv(header: true, sep: '\t', quote: '"')
    def samples = sampleRows(rows, bamsRun)
    def report = ["# DRY RUN: full rerun assuming new output directories; commands below are NOT executed.",
                  "# Sample discovery reads existing FASTQ lists. Paths are production paths for log comparison.",
                  "# Commit: ${['git', '-C', projectDir.toString(), 'rev-parse', 'HEAD'].execute().text.trim()}"]
    def add = { label, command -> report.add("\n# ${label}\n${command.stripIndent().trim()}") }
    sheets.each { sheet ->
        def idx = sheetIndex(sheet.name)
        add.call("DEMUX ${idx}", demuxScript(idx, runDir.toString(), sheet.toString(), "${fastqRun}/${idx}"))
        if (alignPanel(idx) == 'SKIP') { report.add("# ALIGN ${idx}: skipped by panel configuration") }
    }
    jobs.each { job -> add.call("ALIGN ${job[0]}", alignScript(job[0], job[1], job[2], job[3])) }
    add.call('INTEROP_PLOT', interopScript(runDir.toString(), bamsRun))
    samples.each { row -> add.call("SAMPLE_QC ${row[1]}", sampleQcScript(row[1], row[2], row[3], row[4], true)) }
    add.call('QCSUM_MERGE', mergeScript(samples.collect { row -> "${row[4]}/${row[1]}.qcsum.txt" }, bamsRun))
    add.call('MULTIQC', multiqcScript(bamsRun, fastqRun))
    def text = report.join('\n') + '\n'
    file("${launchDir}/commands.txt").text = text
    println text
    log.info "Dry-run commands saved to ${launchDir}/commands.txt"
}

def demuxScript(idx, run_dir, samplesheet, outdir) {
    def fqrun = file(outdir).parent.toString()
    def mk = "${fqrun}/.${idx}.bclqc_demux_incomplete"
    """
    ${writeGuard([fqrun, mk, outdir])}
    ${dragenLock()}
    mkdir -p ${quote(fqrun)}
    touch ${quote(mk)}
    dragen --bcl-conversion-only true --bcl-use-hw false --bcl-only-matched-reads true \\
        --bcl-input-directory ${quote(run_dir)} \\
        --sample-sheet ${quote(samplesheet)} \\
        --output-directory ${quote(outdir)}
    rm -f ${quote(mk)}
    """
}

def alignScript(sample, panel, fastq_list, outdir) {
    def bed = panel == 'HEME' ? params.heme_bed : params.pcp_bed
    def ref = panel == 'HEME' ? params.heme_human_ref : params.pcp_human_ref
    def extra = '--enable-duplicate-marking true'
    if (panel == 'HEME') {
        extra = "--umi-enable true --umi-source qname --umi-correction-scheme random --umi-min-supporting-reads 1 --umi-metrics-interval-file ${quote(bed)} --vc-enable-umi-germline true --vc-enable-high-sensitivity-mode true"
    }
    """
    ${writeGuard([outdir, "${outdir}/.bclqc_align_incomplete", params.dragen_tmp])}
    ${dragenLock()}
    mkdir -p ${quote(outdir)}
    touch ${quote("${outdir}/.bclqc_align_incomplete")}
    dragen --intermediate-results-dir ${quote(params.dragen_tmp)} \\
        --enable-map-align true --enable-map-align-output true --output-format CRAM \\
        --generate-sa-tags true --enable-sort true \\
        --soft-read-trimmers polyg,quality --trim-min-quality 2 \\
        --ref-dir ${quote(ref)} \\
        --qc-coverage-tag-1 target_bed --qc-coverage-region-1 ${quote(bed)} \\
        --qc-coverage-reports-1 cov_report --qc-coverage-ignore-overlaps true \\
        --enable-variant-caller true --vc-combine-phased-variants-distance 6 \\
        --vc-emit-ref-confidence GVCF --enable-hla true \\
        --fastq-list ${quote(fastq_list)} --fastq-list-sample-id ${quote(sample)} \\
        --output-directory ${quote(outdir)} --output-file-prefix ${quote(sample)} \\
        ${extra}
    rm -f ${quote("${outdir}/.bclqc_align_incomplete")}
    """
}

def interopScript(run_dir, bams_run) {
    """
    ${writeGuard([bams_run, "${bams_run}/occ_pf_lane_mqc.jpg"])}
    mkdir -p ${quote(bams_run)}
    MPLBACKEND=Agg occ_pf_plot.py ${quote(run_dir)} ${quote(bams_run)} || echo "WARN: occ_pf plot failed (non-fatal)" >&2
    """
}

def sampleQcScript(sample, qcpanel, bam, outdir, expand = false) {
    def configArgs = ['python3', "${projectDir}/bin/qcsum_cfg.py", params.qcsum_config.toString(), qcpanel.toString()]
    def bait = expand ? quote(readCommand(configArgs + ['get', 'bait_intervals'])) : '"$BAIT"'
    def target = expand ? quote(readCommand(configArgs + ['get', 'target_intervals'])) : '"$TARGET"'
    def perlArgs = ['perl', "${projectDir}/qcsum_metrics.pl", sample.toString(), outdir.toString()]
    def qcsum = expand ? readCommand(configArgs + ['print-perl'] + perlArgs.drop(1)) :
        (['qcsum_cfg.py', params.qcsum_config.toString(), qcpanel.toString()] + perlArgs).collect { quote(it) }.join(' ')
    def hsm = "${outdir}/${sample}.hsm.txt"
    def part = "${hsm}.partial"
    """
    ${writeGuard([outdir, hsm, part, "${outdir}/${sample}.qcsum.txt"])}
    mkdir -p ${quote(outdir)}
    if [ ! -s ${quote(hsm)} ]; then
        BAIT=\$(qcsum_cfg.py ${quote(params.qcsum_config)} ${quote(qcpanel)} get bait_intervals)
        TARGET=\$(qcsum_cfg.py ${quote(params.qcsum_config)} ${quote(qcpanel)} get target_intervals)
        rm -f ${quote(part)}
        picard CollectHsMetrics I=${quote(bam)} O=${quote(part)} R=${quote(params.picard_ref)} \\
            BAIT_INTERVALS=${bait} TARGET_INTERVALS=${target}
        mv -f ${quote(part)} ${quote(hsm)}
    fi
    ${qcsum}
    test -s ${quote("${outdir}/${sample}.qcsum.txt")}
    """
}

def mergeScript(txts, bams_run) {
    def args = txts.collect { quote(it) }.join(' ')
    """
    ${writeGuard(["${bams_run}/qcsum_mqc.csv"])}
    qcsum_merge.py ${quote("${bams_run}/qcsum_mqc.csv")} ${args}
    """
}

def multiqcScript(bams_run, fastq_run) {
    // WGS CSVs are intentionally removed before MultiQC; validation removes scratch links only.
    def pi = "${bams_run}/pipeline_info"
    """
    for x in ${quote(pi)} ${quote("${pi}/wgs_csv_deleted.txt")}; do
        if [ -L "\$x" ]; then echo "ERROR: \$x is a symlink; refusing to write the wgs csv deletion log through it (nothing deleted)" >&2; exit 1; fi
    done
    ${writeGuard([bams_run, pi, "${pi}/wgs_csv_deleted.txt", "${bams_run}/multiqc_report.html"])}
    ${params.forbid_prefixes ? "check_write_paths.py ${quote(params.forbid_prefixes)} --tree ${quote(bams_run + '/multiqc_data')}" : ''}
    mkdir -p ${quote(pi)}
    RC=0
    rm -fv ${quote(bams_run)}/*/*.wgs_*.csv > wgs_csv_deleted.lst || RC=\$?
    { echo "# \$(date -Is) MULTIQC deleted \$(wc -l < wgs_csv_deleted.lst) WGS CSV file(s)"; cat wgs_csv_deleted.lst; } > wgs_csv_deleted.log
    cat wgs_csv_deleted.log >> ${quote("${pi}/wgs_csv_deleted.txt")}
    cat wgs_csv_deleted.log
    [ "\$RC" -eq 0 ] || { echo "ERROR: deleting wgs csv failed (rm rc=\$RC)" >&2; exit "\$RC"; }
    multiqc --force --config ${quote(params.multiqc_config)} --outdir ${quote(bams_run)} ${quote(bams_run)} ${quote(fastq_run)} 1>&2
    """
}

process DEMUX {
    tag "${idx}"
    input:
    tuple val(idx), val(run_dir), val(samplesheet), val(outdir)
    output:
    tuple val(idx), val(outdir), emit: dir
    script:
    demuxScript(idx, run_dir, samplesheet, outdir)
}

process ALIGN {
    tag "${sample}"
    input:
    tuple val(sample), val(panel), val(fastq_list), val(outdir)
    output:
    tuple val(sample), val(outdir), emit: done
    script:
    alignScript(sample, panel, fastq_list, outdir)
}

process INTEROP_PLOT {
    input:
    val(gate)
    val(run_dir)
    val(bams_run)
    output:
    val('done'), emit: done
    script:
    interopScript(run_dir, bams_run)
}

process SAMPLE_QC {
    tag "${sample}"
    input:
    tuple val(i), val(sample), val(qcpanel), val(bam), val(outdir)
    output:
    tuple val(i), val("${outdir}/${sample}.qcsum.txt"), emit: txt
    script:
    sampleQcScript(sample, qcpanel, bam, outdir)
}

process QCSUM_MERGE {
    input:
    val(txts)
    val(bams_run)
    output:
    val("${bams_run}/qcsum_mqc.csv"), emit: csv
    script:
    mergeScript(txts, bams_run)
}

process MULTIQC {
    input:
    val(csv)
    val(plot)
    val(bams_run)
    val(fastq_run)
    output:
    val("${bams_run}/multiqc_report.html"), emit: html
    stdout emit: wgs_deleted
    script:
    multiqcScript(bams_run, fastq_run)
}

workflow {
    if (!params.run_dir) {
        error "Missing required parameter: --run_dir"
    }
    def steps = params.steps.toString().tokenize(', ')
    def unknown = steps.findAll { s -> !['demux', 'align', 'qc'].contains(s) }
    if (!steps || unknown) {
        error "Unknown step(s): ${unknown}. Valid: demux, align, qc"
    }
    if (steps.contains('qc') && !params.sampleinfo) {
        error "--sampleinfo is required for the qc step"
    }
    if (!params.dragen_lock_wait.toString().isLong() || (params.dragen_lock_wait.toString() as long) < 0) {
        error "--dragen_lock_wait must be a whole number of seconds >= 0 (got: ${params.dragen_lock_wait})"
    }
    def runDir = file(params.run_dir).toAbsolutePath().normalize()
    if (!runDir.isDirectory()) { error("Run directory not found: ${runDir}") }
    def runName = runDir.name
    def fastqRun = file(params.fastqs_dir).toAbsolutePath().normalize().resolve(runName).toString()
    def bamsRun = file(params.bams_dir).toAbsolutePath().normalize().resolve(runName).toString()
    if (params.dry_run) {
        if (steps.toSet() != ['demux', 'align', 'qc'].toSet()) { error('Dry-run requires --steps demux,align,qc') }
        dryRun(runDir, fastqRun, bamsRun)
    } else {
        def forbid = params.forbid_prefixes ? params.forbid_prefixes.toString().tokenize(',').collect { x -> x.trim() }.findAll { x -> x } : []
        def guarded = steps.contains('demux') ? [['bams', bamsRun], ['fastqs', fastqRun]] : [['bams', bamsRun]]
        guarded.each { g ->
            def x = forbiddenPrefix(g[1], forbid)
            if (x) {
                error "Refusing to write ${g[0]} run dir ${g[1]}: it is (inside) forbidden prefix ${x} (--forbid_prefixes ${params.forbid_prefixes})"
            }
        }
        if (steps.contains('qc')) {
            if (!steps.contains('demux') && !file(fastqRun).isDirectory()) {
                error "FASTQ run directory ${fastqRun} not found; the qc step needs it for MultiQC. " +
                    "Restore it, or pass --fastqs_dir with the parent dir that holds ${runName}."
            }
            ["${bamsRun}/pipeline_info", "${bamsRun}/pipeline_info/wgs_csv_deleted.txt"].each { x ->
                if (java.nio.file.Files.isSymbolicLink(file(x))) {
                    error "${x} is a symlink; refusing to run qc (the wgs csv deletion log must never be written through a symlink)"
                }
            }
        }
        log.info "bcl-qc-nf | run ${runName} | steps ${steps.join(',')} | fastqs ${fastqRun} | bams ${bamsRun}"

        def ch_demux = channel.value('ok')
        def ch_gate = channel.value('ok')
        if (steps.contains('demux')) {
            def sheets = sampleSheets(runDir)
            if (file(fastqRun).exists()) {
                def present = sheets.collect { s -> sheetIndex(s.name) }.findAll { i -> file("${fastqRun}/${i}").exists() || file(demuxMarker(fastqRun, i)).exists() }
                def bad = demuxProblems(present, fastqRun)
                if (bad) {
                    error "Incomplete demux output in ${fastqRun}:\n  ${bad.join('\n  ')}\n${demuxRecovery(fastqRun)}"
                }
                sheets = sheets.findAll { s -> !present.contains(sheetIndex(s.name)) }
                if (!sheets) {
                    error "FASTQ output directory already exists and every sample sheet is demultiplexed, refusing to demux: ${fastqRun}. To align/QC this run use --steps align,qc (or qc)."
                }
                log.warn "FASTQ output directory ${fastqRun} exists; demultiplexing only the missing index dir(s) ${sheets.collect { s -> sheetIndex(s.name) }} (already complete: ${present})"
            }
            DEMUX(channel.fromList(sheets.collect { s -> [sheetIndex(s.name), runDir.toString(), s.toString(), "${fastqRun}/${sheetIndex(s.name)}".toString()] }))
            ch_demux = DEMUX.out.dir.toList().map { x -> 'ok' }
            ch_gate = ch_demux
        } else if (steps.contains('align')) {
            if (!file(fastqRun).isDirectory()) {
                error "No FASTQ run directory ${fastqRun}; run the demux step first"
            }
        }

        if (steps.contains('align')) {
            ALIGN(ch_demux.flatMap { x ->
                checkDemuxComplete(runDir, fastqRun)
                allAlignJobs(fastqDirs(fastqRun), bamsRun)
            })
            ch_gate = ALIGN.out.done.toList().map { x -> 'ok' }
        }

        if (steps.contains('qc')) {
            def ch_rowlist = channel.fromPath(params.sampleinfo, checkIfExists: true)
                .splitCsv(header: true, sep: '\t', quote: '"')
                .toList()
                .combine(ch_gate)
                .map { joined ->
                    def samples = sampleRows(joined.take(joined.size() - 1), bamsRun)
                    checkAlignComplete(samples)
                    samples
                }
            INTEROP_PLOT(ch_rowlist.map { rs -> 'ok' }, runDir.toString(), bamsRun)
            SAMPLE_QC(ch_rowlist.flatMap { rows -> rows })
            def ch_txt = SAMPLE_QC.out.txt.toSortedList { a, b -> a[0] - b[0] }.map { rows -> rows.collect { it[1] } }
            QCSUM_MERGE(ch_txt, bamsRun)
            MULTIQC(QCSUM_MERGE.out.csv, INTEROP_PLOT.out.done, bamsRun, fastqRun)
            MULTIQC.out.wgs_deleted.subscribe { s -> log.info "MULTIQC wgs csv deletion (also in ${bamsRun}/pipeline_info/wgs_csv_deleted.txt):\n${s.trim()}" }
        }
    }
}
