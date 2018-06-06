# AMrefiner

# License

This software was developed at the Leibniz Institute on Aging - Fritz Lipmann
Institute (FLI; http://www.leibniz-fli.de/) under a mixed licensing model. This
means that researchers at academic and non-profit organizations can use it for
free, while for-profit organizations are required to purchase a license. By
downloading the package you agree with conditions of the FLI Software License
Agreement for Academic Non-commercial Research (FLI-LICENSE).

# Introduction

AMrefiner is a software tool to improve output from ALLMAPS. As ALLMAPS orders scaffolds as a whole, this tool was developed to integrate scaffolds, left out by ALLMAPS.
This is done by using a genetic map and additional information about the size of scaffolds and the gaps within them.
The result in an assembly with additional scaffolds, most of them integrated into gaps, according to the genetic map.

Further details:

- [Input](#input)
- [Options](#options)
- [Output](#output)

# Input

You need four input files:

|command line argument	|description						|
|-----------------------|-------------------------------------------------------|
|--am			|the agp file containing the ALLMAPS output		|
|--marker		|a bed file containing genemap marker			|
|--groups		|an agp file listing all scaffolds and their length	|
|--gaps			|a bed file listing the gaps within the scaffolds	|

# Options

|option	|effect							|
|-------|-------------------------------------------------------|
|--log	|write log file with all interim results		|
|--cut	|cut genetic position of markers after three decimals	|

# Output

|file		|description				|default output	|
|---------------|---------------------------------------|---------------|
|output		|agp file with the new assembly		|Yes		|
|statistics	|text file with run statistics		|Yes		|
|log		|tsv file listing all interim results	|No		|


