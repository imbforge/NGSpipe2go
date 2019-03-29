![IMB-logo](resources/IMB_logo.png)

# NGSpipe2go #

A set of NGS data analysis tools and pipelines developed and utilised at the Institute of Molecular Biology gGmbH in Mainz (https://www.imb.de/).

![NGSpipe2go scheme](resources/NGSpipe2go_scheme.png)

## Prerequisites ##
### RNA-seq pipeline ###
A flowchart for the RNA-seq pipeline is given [here](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=NGSpipe2go_RNAseq_pipeline.html#R7V1pk5s4E%2F41rtr9YBeHuD7OmeyV7M7krezul5RAwmaDwQE8R379KwmEOQTGHgx4JsnWxggBUqvV%2Fai71ZqpV%2BundxHcrP4IEfZnioSeZur1TFFkVTHIP7TkOS0xNSktWEYeyirtCu697zgr5NW2HsJxqWIShn7ibcqFThgE2ElKZTCKwsdyNTf0y1%2FdwCWuFdw70K%2BXfvZQsspKZd3a3XiPveUq%2B7TJO2xD5%2BsyCrdB9r0gDHB6Zw35a7I%2BxiuIwsdCkXozU6%2BiMEzSX%2BunK%2BxTsnKKpc%2FdNtzNmxzhIOnywHvpV3f1wUz8h9Wd%2FeH9%2Fe2X7%2B%2FnfAAeoL%2FNaDFTdJ%2B88BJ5D5S6vrcM2A3925Y29TJiZMgvya9l9i97zI6qJaRJ7F28lBEjeea0XyVrn%2FySyT0f2ti%2FzEl6FfphxCqpt%2BwPqRInUfg1HyRCxEs3DJKMo2SdthvGK4yyN7L35Feu5%2FuFl97o9G%2F%2BUn6HjaF6uYwg8ghtK8VOuPYccinRKj6M4%2Bx3PrxS3skHHCX4qXG45JwJyLzC4Ron0TOpkj1gWdnQZFOKs9FjgT%2F1rGxVYE0VZIUwmxPL%2FNU75iA%2FMv4Q88r3z8n1J%2Bnvf3%2B9vbuJ1jL69O75lznQamOHEZlG2SX27fDxZldQoD0lSRglq3AZBtD%2FPQw32YD8h5PkORs8uE1CUlRgB%2FzkJX%2FTxxdadvVP9jL6%2B%2FqpePHMLwLS1cJD9PIf%2Fj56sXuMXe2eQxdUkJBLh44qHWNaeOv5vDl11vtvu97w7mfs0TDMhGrhNnIysj04H%2F66%2F7z9dXv9%2Fd9Pv3y6%2FSd%2B%2BHeumpngg9ESJy0VreyNlPitnBRhHybeQ1nG9c4WmtrKFrtxH4A35I68UeSMAqP0yRvHM4LMhfJeTuAsMzwntDf8hzY5gTYRsFJNwTRqE1WSy9pErqsTIFl1daIpp1InZyg2yiplJyjEYqNPmaB2lQhgWrpBfu2DPLBu6MwH6qQ0g9qoF6i8LTEI1wT0xjxmQ3xBKshg8yRQE1fK7OJyDb2A0tPbYN8jYpUUXrZrkPSr5eJUQx3bkH66k70l3sCgpgmrBZDygZdUvhMFMMbfFpwWC6L7wofn4qf4SPESqEFHBbqDLVtChg0lzQIWcA3X1SxdUqS5Y0GkYUdzZaQYqoFMQzJ010aupkmSjQEAyIAqtksfWUXYLX1mlSR0MX9BGUm5XXrJamsviPojF97adkPC0OTnh3f3tOnKMiQXth%2FalN1gnGAyire8UzH5fffhgvST%2FNjf4ZRsjfU4T0ABn9SILuQm0XA1MFkVvUyG9eZMWHgJtMnThbmkSAi7XuAlXhgcPKUm3Vs%2FhChmEjDeYIfJwFK%2F1yHa%2BpTXzr%2BLOI6JpvEgldMwQHnn5rTnnkt0FtNGkAh7OtGO6HAPzPeEnW1C6V0YhIaWNGP%2BPqCyrJSgsllHyqohQMpyH0hZqD2BQHtWSEDWCxv6kwjUDVlZ0M6nIInDDpMX5GZMICRXq%2FruTEO5RjJFYKriZZ0xSPa9P0OPcVv2sbmilQbMKL8gdN2YYKXqIOSNfsF61xoT3fLfBRtG75aPPqGr3hW6Tgq56nssGnWJZ4dPVOB5wTIVeXYYIRzNSTG5SMFRNi%2BYQJTIJNDyO0V%2FAr09U1QM6N%2B8xgYilL9badMfNfAiwiliFOkQcRsR2BUvkqeSdaYKHpGlI9nQkG5AGZoSkE1ZIfhRBrZiQqzCueE6kilbLkRAUgwEFQWolqo7kqw5BGBatiZpkkuw4ijgsbGfKZGqt9ugYgMIrA2CWH9y3FVnp2ZeqapWvC4R0SHcTRrlMGYpvFCq1cpZqlBLaearxA7Rc62whgkTVC3xeAEZMPLp6xsGR6D%2FHHtkVMj8lJYpCmNONIlTnzI9HUIpxQSkng7XRFZeBna8SXtB4PwipcJvHx%2Fiz59IV3767eP886efFwW6e4KxqLeyXiLEu%2FVqtLBOhUPo8o7M%2FA3pRECQWJxKMzLjglxqVvlwqh3q%2BP430jBaKJoytDyd9g24luItz4H%2BRWY1TyhCuOQ2dB%2B79KmQ1HJ9ZoBymc2piBUeV16C74kQoG98jOCmTYN3h3egYkw2gADvCQCfeSqELLe7Jn8gsTK1lIlBsf%2FFOPpo%2F0cjO4iUp26btwO1ypK9EWhpjmkYuqu4rgyApbma6shQsxTLcF3kAnsuE1a2FduyCbiyoQosQwYaVkwXIdM1ZAk7DlYA0sYBWg29TAkkVm4%2FQNZBICsmmMjPTSZHKzqv%2Bj4GRg7GHAdraKK7emx8hCFyyNSm8oS%2BOe2DtMOdzd05FQRZptiurx4uC1CRImaX8nWJBYbsXIQ3PkEqSa9DmL%2BSysITdeowbNaOQZSCcsxAS69grn%2FglgeQ7QNusnQ4ciOXBb1%2BAJjLAyZ%2FgLkuYK6zS3hSdjW52SX8mtGeTdERAVSB6y3b4J5t6UCVgGMZDtCQ5EDFUl0AXZX0AxuWOVc03UCSrcoE6WmOBHVZxciRXMW0TYegQmBiG9hmGYcMBveaupmSqHL3B%2BA7CvDVNCZz1NLuPobRV%2Br6m9Eo9AAumfEswqkoyR2aO7V5T90jmQMui2bjrfRDh3oNi4xKkUYQ0icCjInIbbWw%2FbAEnUXDRjRR9Y9qVL5rYx%2Bq0U9ljlLMHwimO4Lh7uzJIJh95qiNUKWsSfO9INUa0maHSdLyOZsZ9J5euEcYPJlns%2BWC6Y0gjcDYhzz%2B3MWpiAI3avpRKBw21bJVrdZRXU2nPb1pimHUJQNgh2KnPIpl8QCjuENkm21rkokQwUgyUEzDQFCWoW0oLlFbro7wHFvIQsgBWDMBmaHQMA3FIrICSlizXFVBCAELKM7JQVQea7SDUHt7m5KsqVpnUJWh4UPY5DRMUUEzmyh88BBl8GL0kojZW0EHawGj7gfyYIlIpZ4Lvps%2BRvsqjJKSfkqhEANBPx8MggpU3WtaUNpNC9PzCAGgdjQsyEMaFkRrzoqePm0gFf%2FgZCOpQHkYyy84XSSV8iOS6hDABCYGmDoOsiiOcTKI6tIjc4Eo6zWhhEPjLugG9GL076vEVLSTC7rIp%2FHdXTCVBk3JBY5uOxghE6gmthUVmaqMbV1CrjmXXceygeXYElAlaAJLtQG0TNsyZNsEGlZtw1AlCY%2BBqfb2NvNHNlQ71FBVxSFHDg9RFnHH0QGS5iAX6raumKatYMd0TFe1bN2yHBcAPNeA4pqOQRAvdjULyzqBwJJlWw4CCDu6jDQZIdmyRhudls4WBqde6wgj4ti4NuewGY%2FK592aZQ69CH%2FbehFGs2KcXE0wHYU3JwYfy7kONENZaB0B5BGeqY7qShsdLYJD6fjq0KJwVyn%2F0uCJKIYAiqUN5%2FTiT5gQkRmwEkXqfUNpZ%2FehPrL17Zaojr%2BumKglTVGfLki%2Fv22J8EroV1mwcOjnkhM%2Bsqaw%2FU%2B0Ab4XfE3f06vyItXI%2F11y45uzSItznbR3WZ%2FSvnFZLxIuFVmdeaJ%2BZ1fXqt4goVlSg%2BzFUg%2FiWuPrNR4ACuoBoMLkNEdY3Lus9sV7z89lSdnrdO66NASjrQ2Fg7VnK9U0RPooA8UzqU1koPYE57zhgQLWtAaqPfXG6x2o%2Fugvhr2aWrF1q5UNpymnZE%2FtRvEiiuBzodqGVohbvlNJ2aMaRoUp0jf2C63bs%2FJMg2t6gdanzubW2Ubbu9g4itsMrcJtmnR6blPbuW1C6YHG0ShDLbk6jlZ7VMs0ZMNr0Ci6MYxGAdaLNErTQv3%2B08VdeZl%2Bh7NgxM2GRkSeckGeW5NjQrTDV%2BSgQL2%2BVuSNDHh8ThLFzDhkpAW3mHh1gykZjfk9HQ5pF7BUESHUYVcgY%2BavK9ssqACoRhquPYSYmKlQfFdyl5FBnTXkZIzCBKbhjansyOYPaZZ2Sf4jpLui4ka7ppvOtUt5d03%2Bo9Wj5CoMyKuhx0YWE3Z8xHEyq6eR7IEHgCQvyolOlLrVRWgjV8HBTNBxxJuzY%2Bz3pbBMQJlHRTrCg%2Ftxm2y2bINX6CPKMhc1ny65%2Bc1pc%2BGO57DNg4vCvBsek2DZJrVduw%2F05hwZkVtkWJ4H9RauPZ%2Fy2nvsP2D61tkQziBDN8tqinuHRovn5Su8yfC5VGP1CMdbP2mNVxiP2YtbMeO25pzi4yzNSJ08g3z73UcKH4PIc1ZrHCSCvRdTm9MnmL%2BaXAcrw6aH4CuXCc9fP1xOdPKSls25bnJC3yeokEUkpHoKQ2c128XKvnpmNgSJs0%2FGzGK%2F1qiOrcmexXD8Mp%2FDi%2F0HL4xmj2lt92SkWh19R5gukiYq2NKgvOubtJGLuzVqCty7vlnEKy94JjUXK%2BxvcBQv7kaPJPsjRJ5L2Y21XEozaBC5lFCB7GzjhJpQioFj%2BeaFXcRZsiLIaEWGME83lu5naAsvm4p0H2j5IQujqkShaCc7KkF7%2BX5CgXNjMG9NLs2zD%2B9Eudhi2yEE6gXSvvfAhwaHilqOaVS1bkbd4%2B2uXuDiCD8R6eStGVmKJlhE5%2F063TFILVZs9vONTWkIFQuc8j07ghnLntxGW2nxweZaMFQAVQ9ipRozlZ%2BZN6mYqZfvwxlP0HQQG4fLohcIGm0cQdM1HuF4QRMFMHne4LhZwjCDD4McbKckmfAUa5BnhpErvIGHCxTtfAUKkCYoULjrYQqL1emf8gS65kgHvXumX3Yy5DkE5JfUw6CKoHezgVgRAH4wJLe7nhxxxlubhth%2FSe37Ipc%2Fu8N8%2FpkqcDFMttFAeiBr3%2BFqQD8fNVANDVD1CaqBPfu7pyEgmsIKjcP3dpdkzQskh9FRH5jyQPqgSRLYcK3Yj2UR4JC5x9J32t7ys7ccYiNO2ozD57txPvO9Cvt0MMH5rp7DfD%2FdtO16kjNfqY02bd99%2FHKTe6fLs3fntS4bkQtndywxS%2Fs4hCYvNfTw%2BW2ez%2Fyu6nNtinYi3qZJz%2B%2BjdtIeu7fgBdLC6iotBlpH1OKvpFOvI9B2cwcRjJoNSqQGywbOTqGUIoYrhpA75LsRbdnhIsc6H5FTsySNHF0sFjnnufugk1F6aJGjdY14eL0ih2KXyxA9X4UPSlnsRNxyQSYQXLIQ3QdmuqCP0JfgYEld7EOIH%2FpJmpeXNEY5WASlo3yeIkibohVDnVDk1RkYsztvczMHWgV1NGafwxbY8YzZA%2B2JA9ykxTVCxz1xx2uE65sv1zf3%2BFtFHVx7Lg1S4Alo8dMmwnGc4tDdsngQLIrJv4erAXCKPW0DLX7ByMatVmPrZOIu69HkSQSdrxMNu0ytwbTLWSMzk9KDF9NsTnG2yHuDGyV0Valx%2B8k2Oonjf%2FQxEc4QOSD6jS0%2FtRas5tU2QGXkU31d04K1F6n8vCUe1Fp9UU9bzBWupisNbmpXtb5crj97adoJsQCfdCLiO8zwhkPXmu9wEJZOX6zJ8eklIWbb%2F%2Bmx1Qg%2FsY1Cy1InqrH1eUCAxM4v8pj0Z99YJu7oMfa7bK1ZeH0xaSs7mXvfAQyvITtrJT2FVd8CxVIGV9XU6XKzmqPnZs1lyNvNzdq6R%2Bqsz%2FXL9wSxudWOot%2FasXMpScppC2Zp%2FlVIxDcL9SJynC2W813X6d4CyNhHlLM13Y5Q8DKTVfAmXrSKzDd7WNtkGzYb7xS5gfaCKXXNJ1ygGYcrvn2HjqXSk9NkN6H34cotLygcGZYbl%2FLx2QrGrCTZXmzgSiKMhQauvQYsWS4OUt2ClWUQKnCDMNsPZRpeRLTEjU7%2FzgSZhxozDYmYrLz4LZll2SBe5vqp8O1b9mdWTxNErwv10j%2F9sLJaTTKmClCcIkBxPA%2F%2BIPuP2uNKJulhmY1of9h%2FOrScjd6ZnRWUt3sUI%2BsvQZe97fARwQTugYcjLcVZRvrZLqnUT9%2B9zaZ8TN04C%2BoArlOMSMDfA104M0dKagHewOywE5yh%2Ba570AdsfwzXGx%2FLizt5wWi8WH6fdjuVI9p51lb1irWiK2Q7nbGCB0NP2F9EM3OWkj5MSJTZcD2rpA57U14hVTO7HoZ0Or%2FQhLZxvoqcQ3yr3v7gmKGOf2jc40cT39xleXxKkRE7Kw%2B9mVtxUlDE3bnfM3fui9ePXbb7ZTl6aFMPDpOw2heZkwqTMMyKWcKQBSJiyECJ99Kv7uqDmfgPqzv7w%2Fv72y%2Ff38%2Bb4yRq6oXKVqHuqNpFgUhF%2FI6XFDLudEH6upqKsBtVxBE2hKyoT%2FNBAx8IuKXZyyBxxZDnGakzBhAwRh%2FJi4RcUD902aazerZLGtjBcdM1y3a%2FtBMFoYqcNlVnem%2B0e4mjWq8ix4KjoG5GZIQVLSvIOoIu2%2ByY%2FrPZ2RB3ozM6QjwSFh5iXz4Fd8lKx6kJDrfAdWOvU3tR99CrmeXfjhNVSIMz9KFynydDdJlTsGVWvjWvaSsx3qzDcLING9GTeRJNY3bUNH2cFCGUaGMm0PdSa%2F5tHrJAl4ixl4TRcyuvnBeaGdzI1S%2BLKmpdtYvBkH4iFm3OEZ9by39JCD0d8dnrbEC3MdspmHFcnkv3zSw2FW6YGG21Wbe15xt46OeYiY8MBkW1LOUxizwlw4Zm3Atl44SJEcnxYcQCld7M8OXbb%2FYOn3Xw8JHLKKSCbgeUSbdWf4SIWhpv%2Fg8%3D).
#### Programs required ####
- FastQC
- STAR
- Samtools
- Bedtools
- Subread
- Picard
- UCSC utilities
- RSeQC
- DEseq2
- dupRadar (provided by another project from imbforge)

#### Files required ####
- targets.txt (sample names)
- contrasts.txt (pairwise comparisons)
- raw reads (.fastq.gz) or mapped data (.bam)

### ChIP-seq pipeline ###
#### Programs required ####
- FastQC
- Bowtie
- Samtools
- Bedtools
- Picard
- UCSC utilities
- MACS2
- ChIPSeeker
- encodeChIPqc (provided by another project from imbforge)

#### Files required ####
- targets.txt (sample names)
- raw reads (.fastq.gz) or mapped data (.bam)

### DNA-seq pipeline ###
#### Programs required ####
- FastQC
- BWA
- Samtools
- Bedtools
- Picard
- dupRadar (provided by another project from imbforge)
- GATK

#### Files required ####
- raw reads (.fastq.gz) or mapped data (.bam)

GATK requires chromosomes in bam files to be karyotypically ordered. Best you use an ordered genome fasta file as reference for the pipeline (assigned in *essential.vars.groovy*, see below).

## NGSpipe2go preparations ##

### Put NGSpipe2go into the project dir ###
NGS projects should be run in a consistant and reproducible way, hence NGSpipe2go asks you to copy all tools into the project folder, which will ensure that you always use the same program versions at a later time point. This can be done either from a local NGSpipe2go copy, a version from the GitHub releases (https://github.com/imbforge/NGSpipe2go/releases) or using the most recent development version from the GitHub repository

    git clone https://github.com/imbforge/NGSpipe2go.git <project_dir>/NGSpipe2go

### Choose one of the pipelines ###

Select a pipeline to run and make symlinks in the main project dir, e.g. for RNA-seq project

    ln -s NGSpipe2go/pipelines/RNAseq/* .
    ln -s NGSpipe2go/modules/RNAseq/essential.vars.groovy .
    ln -s NGSpipe2go/modules/RNAseq/tool.locations.groovy .

or for single-read (SR) ChIP-seq project

    ln -s NGSpipe2go/pipelines/ChIPseq/* .
    ln -s NGSpipe2go/modules/ChIPseq/essential.vars.groovy .
    ln -s NGSpipe2go/modules/ChIPseq/tool.locations.groovy .
    
or for paired-end (PE) ChIP-seq project

    ln -s NGSpipe2go/pipelines/ChIPseq_pe/* .
    ln -s NGSpipe2go/modules/ChIPseq/essential.vars.groovy .
    ln -s NGSpipe2go/modules/ChIPseq/tool.locations.groovy .

### Customise NGSpipe2go to your needs ###

Adjust the project-specific information in the following files:

- *essential.vars.groovy* specifies the main project variables like project dir and reference genome
- *xxx.pipeline.groovy* describes the pipeline steps and the location of the respective modules
- *targets.txt* and *contrasts.txt* contain the sample names and the differential group comparisons
- *tool.location.groovy* and *bpipe.config* specify the paths and resource allocation for the tools

Additional software parameters can be customised in the *xxx.vars.groovy* files accompanying each bpipe module.

## Run a pipeline ##

Copy the input FastQ files into the <project_dir>/rawdata folder.

Using GNU Screen (for persistence) load the bpipe module customised for the Slurm job manager, e.g.

    screen
    module load bpipe/0.9.9.3.slurm

Start running the pipeline of choice, e.g.

    bpipe run rnaseq.pipeline.groovy rawdata/*.fastq.gz

or

    bpipe run chipseq.pipeline.groovy rawdata/*.fastq.gz    

or

    bpipe run chipseq_pe.pipeline.groovy rawdata/*.fastq.gz

## Compile a project report ##

The final result of the provided pipelines will be saved in the ./reports folder.
The Rmd file can be edited or customised using a text editor and then converted into HTML report using knitr
    
    R usage:
    rmarkdown::render("DEreport.Rmd")
    or
    rmarkdown::render("ChIPreport.Rmd")
