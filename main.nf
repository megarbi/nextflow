#!/usr/bin/env nextflow

nextflow.enable.dsl=2



workflow {

    //  Input data is received through channels
    pppp(params.stringa)

    pppp.out.view()
    
}



process pppp{

    input:
    val stringa

    output:
    stdout

    script:
    """
    echo ${stringa}
    gatk --version
    """
}

