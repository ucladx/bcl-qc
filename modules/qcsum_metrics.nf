// modules/qcsum_metrics.nf - Parse Picard HsMetrics and apply QC thresholds

process QCSUM_METRICS {
    tag "${sample_id}"
    label 'low_memory'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(panel), path(hsmetrics)

    output:
    tuple val(sample_id), path("${sample_id}.qcsum.txt"), emit: qcsum

    script:
    """
    #!/usr/bin/env python3
    import yaml
    import sys

    # Load qcsum config
    with open('${projectDir}/config/qcsum_config.yaml') as f:
        config = yaml.safe_load(f)

    # Get panel aliases from config (centralized mapping)
    panel_map = config.get('panel_aliases', {})

    panel_name = '${panel}'
    panel_key = panel_map.get(panel_name, panel_name)

    # Handle heme subpanels
    if panel_key.startswith('heme_') and panel_key != 'heme_comp':
        panel_config = config.get('heme_comp', {}).copy()
        subpanel_config = config.get(panel_key, {})
        if 'target_intervals' in subpanel_config:
            panel_config['target_intervals'] = subpanel_config['target_intervals']
    else:
        panel_config = config.get(panel_key, {})

    # Extract thresholds
    pipeline_version = panel_config.get('pipeline_version', 'unknown')
    platform = panel_config.get('platform', 'NovaSeq6000')
    pass_min_align_pct = float(panel_config.get('pass_min_align_pct', 97))
    fail_min_align_pct = float(panel_config.get('fail_min_align_pct', 97))
    covered = int(panel_config.get('covered', 250))
    pass_min_roi_pct = float(panel_config.get('pass_min_roi_pct', 90))
    fail_min_roi_pct = float(panel_config.get('fail_min_roi_pct', 90))
    pass_min_avgcov = float(panel_config.get('pass_min_avgcov', 400))
    fail_min_avgcov = float(panel_config.get('fail_min_avgcov', 400))
    pass_min_reads = int(panel_config.get('pass_min_reads', 66000000))
    fail_min_reads = int(panel_config.get('fail_min_reads', 66000000))
    capture = panel_config.get('capture', 'unknown')
    capture_version = panel_config.get('capture_version', 'v1.0')

    # Parse Picard HsMetrics (line 7 contains the data)
    with open('${hsmetrics}') as f:
        lines = f.readlines()

    # Find the data line (after ## METRICS CLASS header)
    data_line = None
    for i, line in enumerate(lines):
        if line.startswith('BAIT_SET'):
            data_line = lines[i + 1].strip().split('\\t')
            break

    if not data_line:
        print("ERROR: Could not parse HsMetrics file", file=sys.stderr)
        sys.exit(1)

    # Extract metrics (Picard v2+ column indices)
    TOTAL_READS = int(data_line[22])
    PCT_PF_UQ_READS_ALIGNED = float(data_line[32]) * 100
    PCT_SELECTED_BASES = float(data_line[6]) * 100
    PCT_ON_BAIT = 100 - (float(data_line[7]) * 100)
    MEAN_BAIT_COVERAGE = float(data_line[9])
    MEAN_TARGET_COVERAGE = float(data_line[33])
    MEDIAN_TARGET_COVERAGE = float(data_line[34])
    MAX_TARGET_COVERAGE = float(data_line[35])
    FOLD_80_BASE_PENALTY = float(data_line[44])
    PCT_TARGET_BASES_1X = float(data_line[45]) * 100
    PCT_TARGET_BASES_20X = float(data_line[48]) * 100
    PCT_TARGET_BASES_100X = float(data_line[52]) * 100
    PCT_TARGET_BASES_250X = float(data_line[53]) * 100
    PCT_TARGET_BASES_500X = float(data_line[54]) * 100

    # Determine alignment QC
    alignqc = 'PASS'
    if PCT_PF_UQ_READS_ALIGNED < fail_min_align_pct:
        alignqc = 'FAIL'
    if TOTAL_READS < fail_min_reads:
        alignqc = 'FAIL'

    # Determine coverage QC
    cov_th_map = {20: PCT_TARGET_BASES_20X, 100: PCT_TARGET_BASES_100X,
                  250: PCT_TARGET_BASES_250X, 500: PCT_TARGET_BASES_500X}
    cov_th = cov_th_map.get(covered, 0)

    covqc = 'PASS'
    if cov_th < fail_min_roi_pct:
        covqc = 'FAIL'
    if MEAN_TARGET_COVERAGE < fail_min_avgcov:
        covqc = 'FAIL'

    # Write qcsum output
    with open('${sample_id}.qcsum.txt', 'w') as out:
        # Header
        out.write('Sample,Sequencing_Platform,Pipeline_version,Alignment_QC,Coverage_QC,')
        out.write('Total_Reads,%Reads_Aligned,Capture,Avg_Capture_Coverage,')
        out.write('%On/Near_Bait_Bases,%On_Bait_Bases,FOLD_80_BASE_PENALTY,Avg_ROI_Coverage,')
        out.write('MEDIAN_ROI_COVERAGE,MAX_ROI_COVERAGE,')
        out.write('%ROI_1x,%ROI_20x,%ROI_100x,%ROI_250x,%ROI_500x\\n')

        # Data
        out.write(f'${sample_id},{platform},{pipeline_version},{alignqc},{covqc},')
        out.write(f'{TOTAL_READS},{PCT_PF_UQ_READS_ALIGNED:.2f},{capture},{MEAN_BAIT_COVERAGE:.2f},')
        out.write(f'{PCT_SELECTED_BASES:.2f},{PCT_ON_BAIT:.2f},{FOLD_80_BASE_PENALTY:.2f},{MEAN_TARGET_COVERAGE:.2f},')
        out.write(f'{MEDIAN_TARGET_COVERAGE:.2f},{MAX_TARGET_COVERAGE:.2f},')
        out.write(f'{PCT_TARGET_BASES_1X:.2f},{PCT_TARGET_BASES_20X:.2f},{PCT_TARGET_BASES_100X:.2f},{PCT_TARGET_BASES_250X:.2f},{PCT_TARGET_BASES_500X:.2f}\\n')
    """
}
