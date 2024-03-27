process PLOT_VF_HEATMAP {
    tag "VF_HEATMAP"
    label 'process_single'

    errorStrategy 'ignore'

    conda "${projectDir}/environments/plot.yaml"

    input:
    path(vf_outputs)

    output:
    path("vf_heatmap.html")         ,  emit: main
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p vf_reports
    cp -L ${vf_outputs} vf_reports

    run_notebook.py \
    --notebook ${projectDir}/bin/make_vf_heatmap.ipynb \
    --output vf_heatmap.html \
    vf_reports

    """
}
