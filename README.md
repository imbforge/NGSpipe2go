![IMB-logo](resources/IMB_logo.png)

# NGSpipe2go #

A set of NGS data analysis tools and pipelines developed and utilised at the Institute of Molecular Biology gGmbH in Mainz (https://www.imb.de/).

![NGSpipe2go scheme](resources/NGSpipe2go_scheme.png)

## Prerequisites ##
### RNA-seq pipeline ###
A flowchart for the RNA-seq pipeline is given [here](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=NGSpipe2go_RNAseq_pipeline.html#R7V1pk5s4E%2F41rtr9YBeHuD7OmUnebLI7k61s9ktKIGGzwUAAz4zz619J3CBs7MHAHJvUxohLarWeftRqNTP5Yv34LoTB6g8fYXcmCehxJl%2FOJEmUJY38Q0u2SYmuCEnBMnRQelFRcOf8wmlhdtnGQTiqXBj7vhs7QbXQ8j0PW3GlDIah%2F1C9zPbd6lsDuMSNgjsLus3Srw6KV2mpqBrFiRvsLFfpq%2FWswSa0fixDf%2BOl7%2FN8Dydn1jB7TNrGaAWR%2F1Aqkq9m8kXo%2B3Hya%2F14gV0q1kxiyX3XLWfzKofYi7vccCN8sFef9Ni9X92an27urr%2F%2FuplnHXAP3U0qi5mkuuSB58i5p9J1naXHTqg%2FN7Sq5yETQ35Ifi3Tf9ltZlgvIVViz8pKmTDibSb7Vbx2yS%2BRnHOhid3zXKQXvuuH7CL5mv1HLoni0P%2BRdxIR4rnte3GqUaJK6w2jFUbpE9lz8iPbcd3SQ69U%2Bid%2FaHaG9aF8vgwhcohsa8WWv3YscijQS1wYRenvvHuFvJHlvkm76x6HMX4sFaV99Q77axyHW3JJetYw0q5Jh1SmRg8l%2FVTTslVJNWWQFsJ0TCzzRxfKQX6k%2BsHXlV9f48svwj%2F%2Ffri%2BvQrXIvrybvt%2BDpRG32FEhlF6iF3Tf7gqCkqypyLxw3jlL30Puh99P0g75D8cx9u08%2BAm9klRSR3woxP%2FQ29fKOnRt%2FRh9PflY%2Flgmx14pKmlm%2Bjht%2Bx59KC4jR0V96EzCiTk0KK9SvuYFl47bladpur9t1kHWfNT9Wjt%2BcjfhFYqtnvr0193XzcfNpe%2F%2Fv3y%2Fsv1t%2Bj%2B37msp8AHwyWOd1xopE%2Bkwt%2BpSSF2YezcVzHuKWrBrY74hiFTxRBZEKsYIjZBBAhGE0RU4VQgIu%2FEkAIkBgCS7MxhQFJABx9Imsr1BEjoiAgyGAsRuH2siC%2B9k4%2ByFgPogTwpyyC32gWKtxUFySwBPTGPWBefkQtEEDxyzMSFNDs7X0PHo%2FJ0Auw6BFZJ4fluC5K8tVqcWKhjK9JPc9KnRAH0GpawXgCpHjhx7T2hByP8c5HJYkFsn3%2B%2FLb8q66msBCrQkoFqYcMUkGZCQTGAAWzNthVDFSRhbhkQKdhSbBFJmqwhXRM01TaRrSiCYGIAANKgjM3KS1YhtiuvWcUxncKdUUWSrpdOvNqYC2L%2ByIGzNm2fKDT5%2BendHa26tPTJgen6JlU3GMWY9OJ11qiI%2FL79dEbaSX7sb3AittbrMp2AHD1pCJ2rTbzualGyOnuZjOrNGVg4MTTJ3aWxRGby2HY8J3Z8OsZ%2BM%2BmZ3w8eXJNut%2BtDFDEsjAJsMTQk83naUDp799HGpSr3%2FNuHo4gYHAdSuIYeyhs3p812bGK6mFGCBPPpeDuiwT3o4CO2NjGVd0kHW2rSTv37YMyiVGHMepMwyxqHMItSD4SZa0QBx4jWRECmDQH9SXA1IBMM2viEK2XsQ88Kch8W4IprpxXvLEOxITKJ46fIyjpTkfR9f%2FoO07b0ZXNJqXSYVn2Ab9sRoUz1Tsgr%2FYRprzEmyc1%2Bf6sQ3n0kt0JxC8Y7wExG7cpgJ0Vg1T2OjSbimf4jBTzHWyaQZ%2FohwuGcFJODhCOl44IBIp1mK%2FmZsjOZnp5JMgb0T35FABHKny3tsh8NDsOjK3wyaRG4DQn7ihbxY8VJU%2BeQyFCRqClI1aAIdQGIuigRGikCU9IhluFcsy1BFw0bIiBIGoKSBGRDVi1BVCzCMw1TERTBJpRxFA7Z2s5ESPXTuxhjCxdsdALffmb0q6lO7bpSN614XRGiRbSbVMpiylJ6oNC4Klep0lVSu17Fpo%2B2jcIGIYxRvcTJCkiHkVdfXjE6At1t5JBeIeNTWCYUjK2gCJn0qdLTLhQSTkCuU%2BGaYOW5Z0ZB0grC6heJFP73%2BT76%2BoU05bf%2FfZ5%2F%2FfL7oiR3h9MXzVo2S7hkt3kZLWxK4RC5vCMjPyCN8AgTixI0IyPOy1GzrodTbVDH57%2BSitFC3pCh5cmwb%2BG1lG85FnTPUud5TBnCeeZKd7FN7%2FLJVbbL%2FFA2cz2VucLDyonxHQEB%2BsSHEAa7LHh3egdqPmUNcPgeh%2FDpp2LI4u51qTcmVpWWNDEq9neEw8%2Fmf3RZn6A8Xb15PVSriuytREuxdE1Tbcm2RQAMxVZkS4SKIRmabSMbmHORqLIpmYZJyJUJZWBoIlCwpNsI6bYmCtiysASQMg7RamllIiC%2BcXsjWQeRrIhwIjd3mRxt6Jz68xgZOZhzHGyhie3qsfIhhsgiQ5viCX1y0gah4J3tzTkVBVkm3K6vFi5LVJEyZpvqdUUFhmxciAOXMJW41y7MH0mx8ESNOoyb7eYgUsk4pqSlVzLXP3HLo4f2ETfxiGAAcliy6weQuTxa7o3MdSFznVeGJ%2BVXE9tXhl8y22OLS4RQebaz3EX3TEMFsgAsQ7OAggQLSoZsA2jLpB1YM%2FS5pKgaEkxZJExPsQSoijJGlmBLuqlbhBUCHZvA1Ks8ZDC619bMRES1s2%2BE7yjC17CYbL2WNvfBD3%2FQpb8ZDUH24JI5z0KcQEm%2BoFmYzTu6PJIuwKVBbVktXd%2Biq4ZlRaVMw%2FPpHR7GBHJ3etjePEHPomIjuqj6ZzVyFrK%2Fj9Wop3JHSfobg%2BnOYLLl7MkwmH3uqIBrUtak%2Bo6XWA0hKDhJUj5nI4OeU0vniILH83S0nDG74SURGPuYx59FmA4vcKNhH7ngENTLVo2rjmpqMuzpSZ1Po84ZATuUO%2BVRLIt7GEYdAtxMUxF0hAhHEoGkaxqCoghNTbKJ2bJVhOfYQAZCFsCKDsgIhZquSQbBCihgxbBlCSEEDCBZJydReaxRQaH2tjYRWdtlnUlVyoYPUZPTKEWNzQShf%2B8gquDl6CWesu8kHawGTLqfyI0VIVVaznlvchttKzdKSvgtoUKMBP3etT48qe51LUi7XQvTWxECQO7oWBCHdCzw5pw1O33aQKrshZONpALVbqw%2B4HSRVNJbJFXRX8PzIL42SHJ1%2B6GajZ%2FsEQmzS%2B8qlKLxIFmQqm7H%2Bv7EpMWNB52FIdyWLgvoBdGOCmsat8Jt9apfn9WrUO6kBv2qOi%2BaczK88twhIiKUZU0Ux6LRJ3QPdjkG%2BkUyS9rIBXV10GD3LsxSgbpgA0s1LYyQDmQdm5KMdFnEpiogW5%2BLtmWYwLBMAcgC1IEhmwAaumlooqkDBcumpsmCgMdglntbm67Ktlx2qLuuzsaO7B5iMqOOvQMExUI2VE1V0nVTwpZu6bZsmKphWDYAeK4AydYtjfB%2BbCsGFlUyERAM07AQQNhSRaSICImGMVrv7GhsqXOaVx3hSh2b3RcBgAxrEnLEmjXLNivkSjhLFzpD%2FHPjhJieszZRTEXeACreOnWNc0%2BMQletraJJC6Ujie5jqy7fWCmjM2ZwqBxfHGM%2BvxHFi%2BUvB6xvZHR7jdYfH9Fce%2BmEubL%2FfujNtp3XVNWTUPFDia9WY%2BpZQGYb8a1fD7QnEd82f%2Bk1MXN%2FXTCzQOQkP56RTvm5IVAb08ax8G4G%2BAmmwwcmJ7ZjjdbNdbwfyXN6NbTkMvJ%2Fm5z4aS2S4tx%2B7nXEJIrR6ojhQWHNsqRrhx%2FZ0aWsttgTlo0ifbDQg3FRshl2piGgGbLLzSVzxBpJF%2F8MP2nAc8G0JtZU0OoJwAM6Ag84jQ%2FgUOBRatlFZFHbCTyN6xVt9gTg6aZWe7bpnS5p0fNSqSxf28RUCqj9OmW4KrIn5OxNRVLUMaapItIQKrI7o80rUZFT%2B5wVuda1cm03d0%2Bu4oYh0gYwRJn%2FfNIq1N8CR7%2Bp8bpilDwNjNKUOtMZAKPk3Qo2oaxbvc7WO5uv3mfrT8uDtztKbBpw8OIsiqoNY1GA8SSL0uZTuftydlv1qNziNNI3CGi48Sl9J%2FkiRUSEdrjzBJSk15fzpFUBj0%2F4I%2BmphozkG%2BELr%2BmJJ70xv6PdIRTRgDU8oevAJTGmy8BV9xIdxvUw3rWDEMOcmsSLkttUDPKsJe9p6McwiR1OgCQdP6Rayjn5S0R3QbFHuaQZHZRzsTgmf%2BnlYXzhe%2BTR0GE9i4k6PuAonjVTtfagA0AQF9UsQlLTQcZdfJHBwUrQscfbU8%2FsX6JjabbShTrhiMCAz5s42LDdk76LqMqcNUIFyMmf1q7IgPHiAPLIPT9vhsMQLN0BWtS7S2je08Pdywqb5Rq%2BhmvHpbp2g917TJ86G2KVUVP1qpnKlh1HC5bPZniT0XOhoeohjjZuvDMMZjxlL%2B9zjnZV5xQvZzl8muIZ5N3vPlP66IWOtVpjr3j%2FdMf0CcavIjbJyrC5V7JpzITHr%2BsvJzp4Sc3mmW2yfNclrJCFtSR2CkNrNSsC0V%2B8MmuZ82YIY8Rfghx1DXKyX7k4fpqf0Yv9n7QYzTmzs96TQbUm%2Bw4xnSRNFNiSWM%2FLq6SSi9s1aosHvbxaRCvH25IrFyvsBjiMFrejByj%2B4SPHpurGai4k6WkILsUUkPOgw1LwYr4zqIhbjFeEGa1IF%2Ba5%2FJLNQkdHKb686YfIDdfjxTj2kV2ZD%2FlP36zLWekYbLUmR%2FP0xQWUn9R9e6IAlL2halnEyD737fEeVsezcYgfCQ45ayaWsrMV0RG%2BTjbeUt8UG%2BfZ%2FsAkro1Fs7mOGcJUOU%2Fuja3V%2BGDHLBgqqq0HAKkHsuXfHZxUINvTt7ONByltgbcd8eUJkKKMAyldYwyOh5TQg%2FE2wFE7ljAnDqMRbGsxGdqUP2wDPAyCZBU8HDqU5wsdQJggdGTLCVOYgE7%2F62ig60cFQO%2BrzV05J3f3hjRmJ09su3N%2FnZzh%2BMihP0CvUVZjd1B0%2FXrlaUHRbUYo2ph0d8X3ZL2AF0LAzrAYgtQM2RjGm3AgG5TW73ATpD4fE1QPNZDVCZqgPckYphGX1B8w9TUh1jqilC4OZIragMCEa8l8qCKARYYeS7VrOsuvznKILVhJNQ4f7trzGe51xqmCCQ53%2BTkM9175RNevbPfPJw4cqe8%2Bf7%2FK17erA7ZY9666oUuf1llilpV1CNtdqejhQ1p%2FPkO6bsGVkf1P3OnFs8mmNJ6v6TS7BPb6mvJwjZP5mtAmuIUIhu2%2BJnIFy6yffNg1ZHZ%2FCJAg7w1pzQ7HB%2BP54EPDyTRyMDEXH0bdJnI0PjTH%2Fz7EOB4fTrQZtnMmhr53JzRw6DS7EyjlOPfR9sK%2Fl6oAFGYuBjKU4JLF5t4zHwO9hT4Ee0u6tj4EENFX0mzXpDLSwWCUaMbzBCNliu4GeUIhV8%2FA4915s5s%2B0OSlm8nZvXt6QiZnEh7vzptdp%2BnxBtokPN6XV98vr%2B7wz5opunRsGi%2BRpZTGj0GIoyhhw8VMehBGjMm%2Fh5sgcIqNdAPNl8HILrCdLtnJBHs2Q9jjEFo%2FJhrrmfiMaZPTSqZeqHsnotm%2BonSq%2BQp3Z%2FAcvgMHtKsv3fD2G9B%2B8qTatUz5Wj0Xduek2tkX1LJI2oGSamvVJNnTSKqtTzqp9i1mfMOi89x32PMr31Nt4Pj0EmqznAP0Q%2FQIP7LdSctKI%2BoB%2FXnUgMC%2BSOYw9GfvWMb26IH9RebhNKa%2FnF3YoRxw3ydVjvq%2Bx8RyDddyYhjNTcAs%2FXXdTJ0u07A%2BeqbhHENeb6bhnRuznvWXOvONSGxs7WbRr%2B1DkolIqrkSZkl%2BXkjgm8WDERxnk%2BV8q3eyzYGXzTfZE1Fakibz3yBa7ATLV%2FvhxclWbDbeFyEH2nomNW0ed2pWX5Hpb2o2ag734xzfsxGnZvs%2FhS2mJGIyX5Ls%2BC1skUd%2BhvI%2Fvfe67DWGDwjGcJr%2BJ5bMfVYk%2BfntlxME1W%2FyjTPX8OA6MZ%2FEOt7TOQXzMSfOsQCmXzDBKdHpuid4wPpHcB24WFzcigsm48Xy17TrKR1Rz2ftcKxN5LratNPN47Jo0gm70mmmxMom%2FAlBmQnXs1oqp1flMJcVvetXj07nMp%2FQFrwXkQPGEDvSssyFPVrANUtEcpvmVaksGqeNpKY8cR4yI04G9Q%2FkP7A4SnbXbJB14yxfCq3mwavHRlmmE1891vTanE0TOfAw5PrxjfDBXn3SY%2Fd%2BdWt%2Burm7%2Fv7rJgtu4di8hmmhuMq1G3V3EeCZh494SeliYQeSxzXMg9lqHtKMsCXg4WZvpVahVtRIItuaNJanDFV87EE1CIWpbbOUmooBOIrRRyIZrhY0vy5t0lE9KxK4dfBnD5LxuCE7Xlwgz5ddX2PsTXZPWb9T66yx5D%2FNxFkUMsHyphRkDkGnbGZE%2FwnyjMpC0Tujs8MjKeEhzrdTaJcodRyaWQaQ3tVryMWlnfr9eteWuGJ5hktLqff7Olkl2TEeX9sy0k5hvNp1lMlWbMQFnpPYGL2jjekjXz8Xy8ZMY%2B4kPvxrv%2BYAm2RI0StwazXUU5KbBp1PgdQTqWd7lu7cP%2F4%2BJvK0%2BB9RZx26idiWrVTb8mymr2aKKWXuiNHmmE3v%2BqUfJynDIfXoJUFsiG7vpUSDF5%2FBwjFeUaflOxD2dprx9E77Zrgfvjjvr8y7u0d%2F628%2FfPv89%2FPZbtv33qeaH3ev95grva7OY1nn68WpfMVpun7mKvhIPzBQTajjB1tW6%2BgHvY19gEAgfVSKs0oLM1N9Qsfx2onoDAB62N%2FQ47TqzD2VO5Alrgu5vUsG8yA3xj5Hj9oXSCW1gga8TwwM6T%2FmynPkQKCDIIIc%2FMkSDTOwP8m%2B2F5gg7OFcq8unzLyZ1clO%2FmlNqG7PacbrHDcMlxKnVtoR3IUJ3uw5Mu5sUvyTxtqBqh9Q03h%2BOSzj%2BBW5mTa4Y4%2FckjgKi77o4gUVn%2F4iAL01f8B).
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
