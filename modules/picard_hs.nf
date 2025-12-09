// modules/picard_hs.nf - Picard CollectHsMetrics for targeted sequencing QC

process PICARD_HS {
    tag "${sample_id}"
    label 'medium_memory'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(panel), path(bam)

    output:
    tuple val(sample_id), val(panel), path("${sample_id}.hsm.txt"), emit: hsmetrics

    script:
    """
    #!/usr/bin/env python3
    import yaml
    import subprocess
    import sys

    # Load qcsum config
    with open('${projectDir}/config/qcsum_config.yaml') as f:
        config = yaml.safe_load(f)

    # Get panel aliases from config (centralized mapping)
    panel_map = config.get('panel_aliases', {})

    panel_name = '${panel}'
    panel_key = panel_map.get(panel_name, panel_name)

    # Handle heme subpanels (inherit from heme_comp with target override)
    if panel_key.startswith('heme_') and panel_key != 'heme_comp':
        panel_config = config.get('heme_comp', {}).copy()
        subpanel_config = config.get(panel_key, {})
        if 'target_intervals' in subpanel_config:
            panel_config['target_intervals'] = subpanel_config['target_intervals']
    else:
        panel_config = config.get(panel_key, {})

    bait_intervals = panel_config.get('bait_intervals')
    target_intervals = panel_config.get('target_intervals')

    if not bait_intervals or not target_intervals:
        print(f"ERROR: Could not find intervals for panel {panel_key}", file=sys.stderr)
        sys.exit(1)

    # Run Picard
    cmd = [
        'picard', 'CollectHsMetrics',
        'I=${bam}',
        'O=${sample_id}.hsm.txt',
        'R=${params.picard_ref}',
        f'BAIT_INTERVALS={bait_intervals}',
        f'TARGET_INTERVALS={target_intervals}'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)
    """
}
