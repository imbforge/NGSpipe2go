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
A flowchart for the ChIP-seq pipeline is available for [single-read](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=NGSpipe2go_ChIPseq_pipeline.html#R7R1Zk5s489e4KvtgF%2FfxGE8ymaSys9mZfJXNvqQECFsbDATwHPn1nyQE5hAY2xhwJkmqYoQAqbvVt1oz%2BWrz9C4C4frPwIHeTBKcp5n8ZiZJoizp%2BD%2FS8py2GKqQNqwi5LBOu4Z79BOyxqzbFjkwLnVMgsBLUFhutAPfh3ZSagNRFDyWu7mBV%2F5qCFaw1nBvA6%2Fe%2BgU5yZq1ipq5u3ED0WrNPm1kE7aA%2FX0VBVuffc8PfJje2YDsNWyO8Ro4wWOhSX47k6%2BiIEjSX5unK%2BgRsGYQS5%2B7bribDzmCftLlgRvhg7u%2BNRLvYX1n3d7cX3%2F7eTPPEPAAvC2DxUzSPPzCpYMeCHQ9tPLpDe3Hlgx1GVEw5Jf414r9Tx%2BzomoLHhJ9V9ZKgZE8Z7BfJxsP%2FxLxPQ9Y0FvmIL0KvCCineRr%2Bgd3iZMo%2BJ4jCQNx6QZ%2BwihK1Mi4QbyGDnsjfU9%2B5SLPK7z0rUb%2B5i%2FN7lAcystVBByEYVtptoMNsvGlQLp4II7Z7xy9Qj7JIm4Yuh5glMCnQhPD1TsYbGASPeMu7K5pMtSwJZWR0WOBPjXWti6QpqywRsDWxCp%2F9Y448A9GH3xa%2BfklefNZ%2BOffD9d3b6ON6Hx%2B9%2Fx%2Brqg13EEHLyN2CT0reHy7ayjAnoAkiJJ1sAp84H0MgpAh5D%2BYJM8MeWCbBLipQA7wCSX%2FkMcXKrv6yl5Gfr95Kl48Zxc%2BnmrhIXL5NXsfudg9Rq92zzmvCSPBlzbBKsExabxGXjacOun9t92E2fQZeTRiPg62kc3A9mDf%2Fn3%2FZfth%2B%2Bbnv5%2Fff77%2BGj%2F8O5cNxvhAtIJJS0eTvZEAv5WSIuiBBD2UedwpZMEdjvibh0yVh8iCWOYhYp2JKIJZZyKacC4mIrfykB2TGICRZHcOYyQ71sFnJHXiOoEldOQIsjIWR%2BDiWBV%2FdSQfJS0GoAN5UpJBbpQLhN%2BWCCSTBOTGPKYofo07iEr4xBETV9Ls9XIDkE%2FgiULoIcxWceOyXYKkXy03pxLq2IGUplNui0Pgd2kDBLEoqXzeXmPTB%2F5YZLNbYGkWPDwXh5DBPmuxHNcQHFOyJFM0gCPZOtShYeuObrqC40hzybAcSdNEQbagY5im7bimJhgaMGTRgBqwJNcQVQuUPrKOoFv6zDpJiFH2mpCGdL1CyXprLbBAwxdoY7kBJlH88%2FbdPRm6tArwheUFFiEgECcQ4%2BU6m1SMf1%2Bt33%2FCE8W%2FOkw5hVxzxwzRoJNCwcNGczMPx1WtpHeSOvYtc8oEUAIs%2FHRhjWALHbrIRwkKyNp5ZZE7fxy8aCY9by8ATkx5XBxCm3I5bKeTiRKrPHC2HiG8y58fjGMsSBAgbBj4Tj65OZk2crFIosIGYF5OVt0RE%2B6BBp%2BgvU0IvAs02DCSZpW%2BD01YlEqasFFXhGWdowiLUg%2BKMFc4KhzhWAEBNgdC8hNz1xAbDmTyqQ6UaRVG1pD7phQuuFqlc2cYijWQSRz%2FQ9bWWcVg3%2FsUIEpt7GNzSS0hTC%2B%2FIHDdGKtCVSTkgz7BnG33cpxZec1%2Bfy0psvuU15LqutNkB7BQsmW1XzUdSDP9Xwyjv6z%2FiJMYs0riCyhrnk1M3AqeCEdD%2FirlaVYQOTCa42Z8kapCjPApxyP2sZrfKXqBye2ZJEOF%2FM17hMBx8ndLbQKiu46YwjxeJE8lV0pVL1RtQ9c1V3JdUVFM1VVlWwSqKZm66zquYs0Nx7EsTYeu4ErAEVxTsAwXCrYqy5rm6oYkAyjaNhxJL2yYZgqh8s1G%2Fa9BrdvvZ0pbMz2qTjbNNFGVkXBTgp%2BNKRYPyaZEUXihUOuVk06hl9RMP4kVOM%2B1xtpEEydref%2BpFRy7jvVHsxZUfBdeDRsyf%2BRBHysfhVehTq8njfXhtk6g%2BJ2eJhGDTUhxvvfd55rV7af%2Ffe5xUn64TRqRM%2FjEesYYm9zoSCNCIDx2FuxhNnghcAsz6mPgpJHHHCiPpQyuqIq3K0lSQXozraqgMBGdEtnAe80c%2FwnRgpZZGMCDLvlMgHu5HvWhudRttuxX5VfKvu88WFbSXzkK7DEqP74sKB4HaJt5cPi3ttlF2%2BzsCJ2UH1RsdoT%2BQuoodbFgjc930WrWoo9apqbIgmKbuq2ojmADyZRdBbgyHjbUTWOuybYuuZaqQ8E1TFNxNdl0XNmVsYlumLppSoqu2k5JBxpQH22aZwqjyt3fGukRGmlNMlKnJZnuYxB9J%2F4v4uwCPlgRNBFHH2UccR2I98RHwLxQLGKbjdILbOI6K9IpEbh%2BQJ7wIcQMdtEmd4%2FWWjo%2BeLDu8Xtgpys9E9dh8nj%2BPh1GO5fXUjJ%2B6yvd9ZXMpzsZfWWfdyzkipQNHj7yU6khhDsNJG2f05VB7mmFe5jAkzlbLa%2Bp3PDTMERTyk%2FW8GkXq%2BJFL6xOzCGstq1rvY6aarrsyU2DrzQtqbp1qOqUh3IWDyCKu8R6LVUwHEcWZVGRDF13gCgCC2tNWGy5mgPngibZsuoqmuIIkiiqJtQEV3NtQZc1AQpQ1GXdMkXx7DpUHnAraFB7p5vCrKlbXadiWu7hsd4anZyHKirqTBgFD8ghFF6O4QHvOUZxA%2BG3KiB0MBTSt%2FjBErxKsOEMoTVoKLxKlSKqDv3RdTQ88O51KUjtLoVDBO%2FjGiXwHquy5NnHCITnEca62tWhoA7oUODZmoPGFbMPTjawqJTQqJRf0EtgcXkjilern0jZ3MjO3bWz%2BfjkzPVfXXMqJc0OnSHX2TOk9a1q0UfxXMFzoUNIiC5upkJdrmT360qR5vb2V3ShQqPpCDpSbJMeeI3F9t9XVGZgOMlPrzFSfmwx703I5OwA4zsghE590xF4pHCi6ShkbFht%2Bz5jyRs9Kg64G0kWwTd%2B2Iu0OVcG9oqVlDAaxQqPFVZEDfOJfKRXb2StQcDQFHL2YqEHaaNKZplCFKXGJrkbQI6w%2FbpIG36m76XwtDqvKXGrExiP0pHxKL3beEcxHrWyJUAW9VbGU%2Buv6rMTGE83shLNVro6306jyyKpbJPlxEhK0U6STd1IZE%2Fg7DeJMK5jTpNEpCFIpH0bygshkZ4xX0etXEGtXEnVTGmUPbXD%2BcmCSx9AEGXW9KRJqD9Pd7%2F7WbvyKHkaPEpXq5pOvzzK0rfOu%2BD2nxBu5n%2FfmYlpv9Pm4sVslevVWu8svkbb2MjHVt2X9ovZO%2BeWFoouc82P%2FBUp7dSkRe1Fsll%2BkaYMJHYOtJe0alWOffM4TTNqcuxgOZUg%2BC2GZd%2FOHWS5FGFIEjrO6cUpJNDQsYiHu3KUApr6cuU0roTj9xZJBiPFkTw1fOAxV00hLkAwMr8nKBEKQddXMYYQzS0inr0%2FagyPBHYLkGVx3bL%2Fi7CYav7EBjkO5YcVJOxa7hhk5FlDNYUoSECatJEyObZ28bDUJf6HoXlF%2BKKKh3mFr8XdNf5HukfJVeDjVwNEkQ0xlT7COJnVC0D0QBaKIC7Ke5ikugePGzCSlYPpoiMRaDUi6B74pJv8WPhTOCLS%2F9c2SbO0SXUmQjKva7F%2FfPOH3RbsHy%2B0nzIrfDfIp4EoY8PTiUrj7hIJPT3PqEiwWQWTa7BBHqG1G%2Bg9QPLW2RBxUV0zSmSuZBWIRstSykzQydC5UCP1CMZbL4mnSeyvMJHhudGQvxBC8H1H7sOmTrxBrgujPEWCDcVDccJJ9pzacjvD0lLFumrBXVrG2ZaWMfml5QWria4rPLJ5JjbswPOwDgedXIRAYK9nWULOCyBmfUg5wQ9fjhq%2FnGxZu0Zc73XsZBjdX8Ou9zSM04oxNJewm4piHEFiv0yUsaUZiMSgTIe5uNs4TXmKpNciXiP%2FGfddrKEXwihe3I2elfln4CCXkBwdO9V3KG9KCFO2t3FCHB%2FcLM14llVgSdZYr1tjNMYplWapk22Zk1Ph8ANZByI33a%2BnfZgd2f7pmxg4kZLBoj05R2cf3rHzswYUz5TAsjfVLcs42efZPd45ugSbzklvu%2F3zs3N6S9OcN%2FyxI1LelKFS3nrgDtUst7yS%2BKSy3NrTkabNL5qycjsyjxP4hToOv%2BiagHA8vwjXwMfaAPVGlLgG5iJY%2FcDvQ3acWqZgky4UwUI%2B3TOLGcoa6w3Y8ItQnNB%2BgwRdCmM%2BnKGol8tQFGFchsKPIytj2p1Hc49euYPW0VLM4rqj7Y3EVhxklRTyhf42bUtjZKQ79FfUTsAwJ1YAsRzAA1ntbgRWG5jWUSSIOb%2FSkI738EWuXc4ir0Zc8%2BD%2FlLQGqV1r%2BNVyx%2FqyMvSOfMEY6hSExpQKsLkiljc9WqbIHGzMFBJahSqINniNxNTHa6HVIykJIiRY%2Bn8fYscMHqDNBng4N9AvhxtURf4kuYF8CdygVwHf9TgTeWwBT3Tib8D3SXYJqYtdXsyYwgFW1Wm97NTJ94Bi4h34mfYeSH0vj%2FHw9WxcznquSnd1ZJ8Adzfv6S6BX97%2BP09a9177X62mlPZv%2F1t2mUl8urrDPZZBkmBDANokR%2BEqgK6LbATZRwdgEdYRrkLzcthCzbIfOc2SyxZGPd6oR8N%2BH6M4ni2cadNi5x3zfSeI19iPedK%2BpCaeg0KMpNTLUGI9oUcL1qU1i7MewzCc3ZAO5juqMLtYvqNqUzQvJpTJMv1z1bJtRfudj8ZAtkk36TLqZrFLK77XGcn9G6BHCRLFqAqs9p1G1f7qaZUZmgTPBtixVJE5aRquDTxvuF1DdByHK7jn2DM0kN2rjOzHanW7TiZTrp7%2FmzlYp5goZ6HVl5IXmOX%2BMp9SzHxKLzC1nee1PVs2MDcqezHFjKYgXY3O0nUsDap12Eexr7PXvr2DdLcLOWpJeAf9oHR6To2LTa%2Fy7fKvL5%2FfvyU4wIviiYa9VqVpVHOBs7sCqe5GauNiWb9NATB0NjAdPwpIAdRCLIAsDvs7CfGl57J03nY0SjlZlrNMStuT0ccUE2Tx7amke1Q114kVZy3rThpnb4lomnVpcr4DHo36jusa4M5biTVndi%2B3Emvr5pNLPI8llwKFzRZ0dbWruy%2FtCJEUJOXN2rM0mRsgP03EdJFPd5XmW3HTDG9e6nca%2BC3s%2FcCGahgvWtnliz1yY7IDm413FshAm2ukutTj2lDVWEh%2FW2tGrXJ9nB96dpARNfCOSjE7j%2FjCDK183KM4irIDIdt3VIJHBxCrY4qOIlruerarMvLqJwrD8hkMv4onqHKuQlcedr4SB%2BaIlNvNxUkKdpV2lU6Icnc75F5mjQ5ZNRYdTwY5X2GDUTeY%2FIKFDUyxoxTO%2FA6jZbPSnfV3rFhAKWzHJknsmtRbRCMPeFF%2Fd4JHf5aVGJgNE8%2FLKgCQcR4c1jOLQJ14WE83Kjq6LnL4w5CBvRvhg7u%2BNRLvYX1n3d7cX3%2F7eTNvjuvVZAthrFzBUXUPKDz58BGu6N6pXBCkr6vJB6tRPrAShAXOwy0XSMRCpalWtbCxSiGPGMoMsgfSEIVMKuTbn%2BuEobR44E4RHFwqqJ8eRQ%2B0ne3KEnXwYA5SdbMGO15aFs97qZwLdqeElrSq2ljwl2Xg3DVSwPKcZel2Y9%2BKyX%2FhrqjnDjujq4dH6oQDH7xaoy5R6rg0lcMPe%2BtGXkOGE1rp%2B%2BVGE7hgucBgAvN2Xqde8Zb1%2BNLCBr%2F95hc1sBEd%2BmeRMUZHGdNHgWguLxuzbi5KfbbXQcUDNslslxfg16qRpyTXBTpfBdLORJ7NtWfzDJb3CYanXcz%2B3CW3UIRuY3qwMaO2vD7fizExpcwdMZqNWXevvwmStBAuIC69NG3JwTCjVRR48Xgafn9BSMtTw%2FcizTwdadzk2IupWDRwciy%2FvtOsm7N4qJBtk6%2F4JBnZb%2BrrEmFzziUFU9KiaEISBF5rEv%2B0FILlbBemw7aMX%2F0qSW5NUFL5KpnkIk%2FX3CVm7qwetiizFhUYgqvYmmVDxzEU2YCWJDuGLEJLExzXmNuaA1xDdTRFF4Fo2LIu6rZlAsu1dV01dMw1VKBbTukj6wi6pc%2Bcyee%2Fd7rMPG3oliEZdCIFHhaGTSc%2Bkh6wnIg7kgNGpu24QLM0LNYtCdqGbbiyaWmmabuKAueaBg1RM3Vbt12oSSK0Jai6iuGaluVi8rAE27VNRxqPHFpmW6CGeq9GYmhAfUuO8y6Bmqs2Vl9QSAMkY0sle8zyrdNMwWIGdnpWOvyxRRHVcfJS0lUOd1QWdnudwYI4qUfJJpC%2FVkn9UMW6bsqvDS3IBwvBLpExLhQljtVxLsdrmzrxgh2vfLAIl6KT9qlbdi6rNdTW9Sblkh1wQsNQH8mRLJW6WuEznUZMT%2FShR7ak1TNz%2Fw9rzNxAZFhnykrYoJh4l4EPgy25ZkOnoc88OUHipic04%2BgyshMUSamcUsc7lmX8irnamGtdWBxkgeKLTzBC9GyImlU6SNULPgw5mcRtLGYqh%2B4eEvjcRt7zkmytJhPkrpkCenf0kV5l9f%2FezM0DYX%2FAejPVynpTZY2XD8TZq4ZtuWNWHOZaSVHyYjis%2Fwwcwqff%2Fh8%3D) and [paired-end](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=NGSpipe2go_ChIPseq_pe_pipeline.html#R7V1pk5s4E%2F41rsr7wS7u4%2BN4JpOjsrPZmWxld79sCRA2OxgI4Dny61%2B1EJhDYHxhZiZJqmKEAKnV6n7U3WpN5MvV04cYRcvfQgf7E0lwniby1USSRFnSyX9Q8pyVGKqQFSxiz2GVNgV33k%2FMCvNqa8%2FBSaViGoZ%2B6kXVQjsMAmynlTIUx%2BFjtZob%2BtWvRmiBGwV3NvKbpd89J12yUlEzNzc%2BYm%2BxZJ828g5byL5fxOE6YN8LwgBnd1Yofw3rY7JETvhYKpLfT%2BTLOAzT7Nfq6RL7QNacYtlz1y13iybHOEj7PPBR%2BOwub4zUf1jeWjcf767%2F%2Fflxmg%2FAA%2FLXjBYTSfPJC%2BeO9wDU9b1FQG9oP9bQ1HlMyVBckl8L9j99zIrrJaRJ9F15KSVG%2BpzTfpmufPJLJPd8ZGF%2FXpD0MvTDmFaSr%2BkfUiVJ4%2FC%2BGCRCxLkbBinjKFGDdqNkiR32Rvqe4sr1fL%2F00vca%2FC1emt%2BhYyjPFzFyPELbWrEdrjybXApQxUdJwn4XwysUnSyPDRuuBxyn%2BKlUxMbqAw5XOI2fSRV21zTZ0LAplbPRY4k%2FNVa2LLGmrLBCxObEonj1hjnID8YffF75%2BT29%2Bib89c%2Fn69v38Up0vn14%2FjRV1MbYYYdMI3aJfSt8fL8pKNEeSBLG6TJchAHyv4RhxAbkP5ymz2zw0DoNSVGJHfCTl%2F4Fj89UdvU3exn8vnoqXzznFwHpaukhuPw7fx9cbB6jV5vnnAsQJOTShlGFMYbCa8%2FPm9Nkvf%2FWqyjvPmOP1pFPwnVsM7I92Dd%2F3H1ff15f%2Ffzn26dv138nD%2F9MZYMJPhQvcNpR0WRvBOJ3clKMfZR6D1UZdwhbcJsj%2FpIhY5UhsiBWZYjYFCKKYDaFiGqcSojInTJkIyQGECT5nd0EyUZ08AVJk7kOEAk9JYKsnEsicMdYFV%2F7IO%2BlLQbgA3lUmkFu1QsgbysMkmsCuDFN6BBfkAqiEj1x1MSlNLmYr5AXAD29CPseEaukcN6tQbKvVoszDbVvQ9hLkggF9TIEA%2BaltdfaS7KkwT%2F%2BjfAsb%2FiMKKrw4bn89pyseYmouo7uKJrlGrrpCKKha0hCqopUjAXVNKaShnXTkFVXlbGiyQpCsu6oumHJsmToouIYiqHLkl35yDLGbuUzyzSF9dYFjLp0vfDS5dqaEV1FLryV5YaE%2B8jPmw930HRpEZILyw8t4A2UpJiQ%2FDrvVEJ%2BXy4%2Ffc36Si76dTwjYWfdfDARZ4B547AdSxyREQ7j6ymdul6KLPJ0ibPJuhq7XuClXggc%2F86CO%2F%2FbmdVH3W8%2FRE5CJVMSkRU%2ByCayuoaOwlo6dNY%2B8NTL7x9OEiL%2BPQTCEwVO0bkpdNtziSKhKgIRCQwTao8OH4EHn7C9ToHeJR5saUk7ED8GfhWlCn41mvBV1jnwVZSOAF%2B5Kk3hqLQaCQiIj%2BAnEZwRgfvQ%2BQy55FjAyAsKi5LCJVenTu1NQ7FBMoljNcjLegMD9r2voUe5jX1sKqmVAdOrLwhdNyEApj4IRaMPWIR22yZODDnz339X4Oc2yFkBnBv8OcC6Ip9W2wHlQHjyzwTHv1v%2FgWmXiEpYwVfxYpsQt8InkGhesMhkmhXGDo6npJhcZCiHMT6VeALhcrW4U7bdwm2ymscK%2FC1qRMhxindLx0GAGc2TWfpUMYDUIZ9qG7quuZLriopiArCzRaSakqm7ruMq1tRwHMvSdOwKroQcwTUFy3CxYKuyrGmubkgywqJt43NAPvKrpZsZhao3T4XochzVZJt2nqjrSLyq0M8mHEuaZFOmKL1QaNQqWKdUS2rnn9QKnedGYaOjqZOXfPraSY5NxeajeYlXfheZDSvov%2BfjgICP0qu8Xq%2BHwmZzOztQ%2Fs6ROpGgVUTHfOu7T9Wrm69%2Ffjtip4JonbYOzuAdO%2FKIsc6dfdBACUT79oI9zBovhG6pR8doOBTyhAOVsVTAlaF4N0iSStqboaoSYAJM6dnIv2Dm%2BhRQ0Dw33vvYhc%2BEpJbrU8uXS41d8%2BNCfqVqsS5cXBX8ygGw%2B0B%2BclkCHjugzcKl%2Bwtt9kGbvc2Xo7Jeiu3my1cER6mJhSC%2BwPUWkw48apmaIguKbeq2ojqCjSRTdhXkyqTZYHucarKtS66l6lhwDdNUXE02HVd2ZbJEN0zdNCVFV22ngoEGxKNt%2FcxoVLv7C5HugUgbmpEaLaG7j2F8D%2FYvMHahAC1gmMDQRwVH0iTiHdgImBWK%2BVnzVvqhDaazMp%2BCwg1CeCLAmAjYWZfe3Ru19HxwZ%2Bzxq2GHg56RY5jCC78Nw2inslpKxi%2B80h%2Bv5Dbd0eCVbdaxiKtSVqT5XpBpDSHaIJCsfEpnBtzTSvcIg6dTNlsuqN4IMjdEW6BOXvB146vieS%2BsXsIhqpctG7X26mo27eGmwQdNcwq3doVOhStn9oDipIcb17JUwXAcWZRFRTJ03UGiiCyCmojacjUHTwVNsmXVVTTFESRRVE2sCa7m2oIuawIWsKjLumWK4skxVOFwKyGord3NaNZWrYmpGMrtYAjeWPD45DRcUYMzURw%2BeA5weNWHh%2FznxEtaGL8TgNDGUErfkAcr9KrQhtOETqeh8C4DRRQO%2Fa9va3jk3WpSkLpNCrso3sell%2BI7AmXh2ccYRadRxrra16CgDmhQ4K01B%2FUr5h8crWNRqQyjUn3BURyL5ufYuL%2B%2Ftu%2F9b%2FMpuvogYvHLG3AsVmJd23BUDePsGcbW2xCkDQukrone%2B%2BOSCl3yZfnpgvT7x5oIrxS%2BapMVeBwCp1DjboweaVNoPAc0gOCe%2BwmLfjii5iXVINqC3Phhz7LiQptulcsZqVvlMk%2BW1GQ1Myp8oVdXstYioWnkNHuxcARxrUpmZaJLitKQM9x9D3ssnvqIa36A60uRCc3pXJntB0xlpedUVo6%2BSKKPXsQxei5ViEBdJO36Q61FwsuiXmaM7fVVvcZIWQv21TVcthLNTr463QaboVmKr2l7clS%2BtXBkHKVowuk5ROrePjEODjmeDee4%2B6v6iizZHAWDKXm0YyGCjstg3EnYXHa8sp0bpwW4vbXiwACXiIXUwxKE4ldA7i1mXpkoAtfQKeFsyRWXNSbCu6Na5QSo9gigtR6nLBn6WUErn3gMtZZsDDAm0zsYFKFkwH0XIS%2FGwBc4APNRTSCAjbhEWGYirq4EYO7WXTErz3GovKiNwabklhFGnrRsp4zDFGX%2Bn2wyMxFKmqXOyT9CzEsQBypp5iW5FjfX5B9Uj9PLMCCvRh4da0zY9BEn6aS5A%2FQIXKEI4qwaDi011zJc25Os7MwWPXlAa%2FBAfxsq3S%2FALKnCHk6D39dpFvAF6RmAZS4abgRy84fd5Tc4n5cgk1VgBCi64VHJRroTV9rdx6h6uMuyzLD5FuZrtPJ84LWP2H%2FA8NZBTKy6ZlShcL7oO5vDU9ZHxudCg9VjnKz9NBkns78jTEb6Rr0HQoTR%2FYbdh%2FXCXHmui%2BPC28Ka4ntJyokbGdt0O8HU0vJptG1qHWMDf2eaihFPLT9cjHRekZZNc7Vhh75PIBwFWpkKwcheTnLf3htgZp2TjuJkeoJvyD2rJXe0eW32XwHnyYe2J7HRuOvucwX75u0ejVRrAuMYw%2FplpIItC2aA9WTWzNntymkLeYBas2TpBc%2Bk7myJ%2FQjHyez27AEev4WO5wLL0bZTvENlUwpC2V4nKVg%2BuAEfySTfzJ0uCa5bkmFMMi7NozC6gjDGIuEHWh2I3MiBI23p6OnU188p9l%2BXU18Zq1N%2Fjla9ffqb%2FXWTU9pAM5c%2B%2BdgeHn1lKI%2F%2BEaZ83Ylf5Ac9kz2UKwTebCay40x7rTnt2926A8A6%2FsffbOD7yQa5w2U6nGwnkpp5kdaB92Nd82plxqsJNdUCHFz7qTctfFwDhG5tfF31hu4u9stac%2BRiv%2B4GkzV91oy2Pbvgz1v5eoXCgOhP7ykhjIH3xhCQdQnLIpr4uywcbDL9U5ptIIxXZKok1ABneYtH2PoppDGy74cI7CQNtFkDdxcK%2BssRCnUsqOY51cckEn7BhIOEQN%2FM1bLOZ5VTCQFwF%2F2LggDCBiB3YlUQkNmBYi%2BhORUz682Dl8AK8WdWexiEUGvj7rLAeDmyoA4QNI4D4OyyQFReuzAYEB6YfSXDwMahTwG5wSykZaGwoKkHPHoXOggVJMGNqem3HBc1tNEoaxK0Z3cJUV6cjVxCqHI16FURRmg5Es9qOnpdEkLt6zAcWkJ4Eek3DhYsLUshISKf5hrJ0s3lNYaBCpsm7SwDVGHycmRAbcWgnDmaln%2BOxYgiB8Z%2FkIXSOzDe4MqWc9mPuzdkvQJBf1r013vHzaHrwqNssFE1scwwW%2BtrykEbctoUzwrZiVTTOVnYo418f7htGrQdu4PNU2zRGGg5qsrNyISzb9sw8m6NJTKpGW%2BZ20zHGJhkeYvvFcMui7Vkpp6EmXreYCixmjtHh4i%2BtPS18yG8%2BSvCq%2Bkft2Zq2h%2B0l7ONfgxZyYzeyvRcCKqz2XuJr5OnLbvFdHcBZMkXPuAgrCQ%2Bb0ix8SUtm%2F%2F%2B%2Fdun99Ri5eAnaptaVLpRj73M7wqQVwTSmhFdv84IMHT0JW2%2FF0LuqpKJHiaHfQ9euyyldu9tHmfJBMZiRCErKbQ%2BoSMBk29LErS9EnGNLK9WzZSvcEz55qBn8xjNDa4Nwp02iVYh7N5uEq3OYP%2BXmEq70AKl4HY6u7rh7lvL%2FpyRpLo5dpLF2SIvyAKvXC%2Bgu%2FiKrY9Z8C0vKjfzx5Zi7clCNUpmneLyzWZLHm3DJudL4zzQZgbOqezcNZR%2Bsh1sZ93KsJ8derLTImrgHWxinvHxhS20inafxVCUn%2BXTvYMNPToIVh1jNBTRRIuTTVaHdz%2B9KKqmzz3PaoOsgjL1SbTjA6wq6Cops2JVIhI2ROetMM7W%2FuwwJHF2K84ojWeLn%2BNup7RHO1%2B0ZbCWIrmvThOFU63k8gxsI7Z5w3aCyrbOEYmyTWDS20ySwd93MHBmgbNGEb7CzAJm3%2FyguSHqbMdy0K3tt2y3fsWPyzoJqjwzH1IlTib1vRM%2BBpN8j%2F9kGAdvvgUf2rmzn9csE3Xkfl7dqC3adPHM%2B5I%2BCp%2Fd5Y2R%2Bg%2FLW%2Bvm4931vz8%2FTtsdvQ3dctgZ71%2FwgkawFoqg5XB5q1U%2FsByAJcnDzdcHaqFW1Egb2JomkMcMVQF5BNYgGKaWZVVqMobCYYxjWLG5XNA8CYIeTjfZ5AXqYdIeJOtlg3YaJ06PZ85WTkW7Q3yNWh02lgyoOTk3hZSwvDUFWUTAms1K4L9ok1RzMzpnh4d7YsKBD1FrcJco9Zyayu4Ht%2FRjryH9S538%2FXbdS1yyvEDvEjN%2FX2duko75%2BNb8SJ3EeLOOlNE27IwenpPoGKOnjjlGhmauLDtn4lovM%2BJfhzUL2CjDn96AXavBnhIvQJkLgbQTsWd78tfCQP4pJfS0y%2BHAG9s5HdB1QpOzMG4rEuS9mSWmlJsjzrbGbJrXr8I0y0SLwKSXxbE5hGY0UwYvQIPGY7yhQSv2CmwdNPPwQeNGS7%2F6bel7Wou5xOprLP51hHNh5Jt7ZDnnQlKclIpvIQ1Dv3NXx7gAwXyyzyHO0MlZEb%2Fb4xBnFRmCq9iaZWPHMRTZwJYkO4YsYksTHNeY2pqDXEN1NEUXkWjYsi7qtmUiy7V1XTV0IjVUpFtO5SNDHeK8tbtsedpSLR9k1IsVznt88wH8QPRE0pMdyGDajos0SyNq3ZKwbdiGK5uWZpq2qyh4qmnYEDVTt3XbxZokYlvCqqsYrmlZLmEPS7Bd23Sk87FDR29L3NCs1coMLUPfEfS%2Biajnwsb6C0pxodC2TLMnLAA%2FCx0th%2BRnx7biH2sWLFPkcq5LuIPOx27XQ61eshEENNZiP3pvChMFeWcl2MczxqWixFl1nMrw2gUn3rDhlU8W4aVg0mNiy97ZzobKZdAGLtkJI9QN9QXORKnlP4ueaTcSeqQOPTOFSMRJKYifFeZmIGjWiaISVl4C1mUU4HAN16zp1PVZBCdI3PCE9jF6GdEJiqTUjonjnYsyZHQCf65r55zrwmynFSi5%2BIpjjx7O0FiVDpIGhU9DTmh5l4gZPrS8s9m99O869p%2FnsNceOsidM6Xh3fBHdpXnabyamjvSfof5Zqq1%2BaYKCi8eiLN5kazlDqcxN9PMiMIFx52N%2FqCzwyvZSQbPJvRiTh3pdcb7sAN3ogTzPRP%2BiS9l5Npsr8r2XVbkqq4zR8gHJ8og3JMPfuUDG2KQDz5Ngj56aDqw%2FLpvOjDZNGoMdpR0YLerq3VUW7r9hmKyZruO8Sp8AKOZs4586qxzJkMcb5FlpI1X5LOdgeNbz4M5%2FcqsAQs5nNuKFOv5wRR10HUZ0DgM0zIDka4tfwsd4I33%2Fwc%3D) sequencing.
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
A flowchart for the DNA-seq pipeline is available [here](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=NGSpipe2go_DNAseq_pipeline.html#R7V1bd5u4Fv41XqvzYC%2Ful8c4aZL2tDnTpD2dzsssgYTNBCMKOIn7648kBOYibOzYhjSZ6UyNuElbW3t%2F%2B6LNSD1fPF3FIJp%2FxhAFI0WCTyP1YqQosqqY5C%2FasspaLF3KGmaxD%2FlF64Y7%2FxfijfllSx%2BipHJhinGQ%2BlG10cVhiNy00gbiGD9WL%2FNwUH1rBGao0XDngqDZ%2Bt2H6Zy3yoa9PnGN%2FNmcv9rKB%2BwA934W42XI3xfiEGVnFiB%2FDB9jMgcQP5aa1Pcj9TzGOM1%2BLZ7OUUDJmlMsu%2B%2By5WzR5RiFaZcbrqWP3vzGSoOH%2Ba1zc313%2Bc%2Bv63E%2BAQ8gWHJajBQjIA%2BcQv%2BBUjfwZyE7Yfxc0q5OY0aG4pD8mvG%2F2W1OXG8hXWLPylsZMdJVTvt5ugjIL5mcC4CDgmlB0nMc4JhdpF6yf8glSRrj%2B2KSCBGnHg5TzlGyQfsNkjmC%2FInsOcWR5wdB6aHvDfpv8dD8DJtDdTqLAfQJbWvNLl74LjmU6CUBSBL%2Bu5heqRhkeW74dD2gOEVPpSY%2BV1cIL1Aar8gl%2FKxt86nhSypno8cSfxq8bV5iTVXjjYCviVnx6DVzkB%2BcP8S88ut7evFV%2Buvvj5e37%2BOFDL9erT6MNb0xdwiSZcQPUeDgx%2FfrhhLtKUlwnM7xDIcg%2BIRxxCfkX5SmKz55YJli0lRiB%2FTkp3%2FR2yc6P%2FrBH0Z%2FXzyVD1b5QUiGWrqJHv7In0cP1rexo%2FV98IwKEnLo0lmlc0wbL%2F0g706T9f5dLqJ8%2BJw9Wmc%2BwcvY5WR7cG%2B%2B3H1fflxe%2FPr764evlz%2BSh7%2FHqsUFH4hnKN1woc2fSIm%2FkZNiFIDUf6jKuOewhbA78psMGaoMUSW5KkPkphDRJLspRHT5WEJE3ShD1kLiBIIkP7ObIFmLDrEgaTLXM0RCR4mgan1JBOEc6%2FLvPsl7aYuKUDkBU6iDUhNqq5KgwrfCLblaoCfGCZvvM3KBrEVPAp1xrozOpgvgh5SefoQCn8hY0jjdrE6yt1abM3W1b0f4Q5IIhPU2QCfMT2uPhSFI0M9J3ukJ0Vj4YVV%2Bck7SvEUGpq24ko6gqWlQcWzLcx2gE%2FvLcBSk6mOIZF0Fuiy7SDE9pLo6BBY0NBUqSEOa7OiGDGzPqrxkHiOv8pp5mlLD64zOuHI589P50pkQpUUO%2FIXjYcJ55OfN1R3tujLD5MAJsEP5AiQpIuS%2BzAeVkN8XN2dknFSFbx1wRrbW6%2FLJA4IJFdG9DhsOOM3P49oxW5h%2BChxyd4lviQmNPD%2F0Ux9Tfn7n0DN%2F7MzIgx53gAFMmNxJImLMU8lDDGk6UGo2Y7gMKNe8%2FPGhJCHC3QdUNIIQFoMb02H7HtEZTAEAIl%2FpktljwAfgwSfkLlNK7xIPtvSkHXMfAqrKSgWqWk2kqpoCpCorB0CqQoWlCRRWjQQEr0f0JxGNEUH2dPAZSMnVvpU3FM4jTUiujRqzMw3lBskUgYMgb%2Bus9vn7%2FsQ%2B4zb%2BsrGiVybMrD4Ae15C4El9EopO7z8vhT%2BtF3SZ%2F%2F5RQZrb0GUFW66h5glMCLkzXBwUWpTb4WKblHfwExV5fjjLhJ6DY4jiMWkmBxnQ4SuDiUSJsJFenCn7celpYtlTvIS04ooIQFg8WzkMAGRKjwCr0PNnm2CfYxMEJ2mubbqaDiUXKLbqacBTSbeRaVtjqJm2ZMmyCW0DWbJpy1B1LVc1Fd3zLBW4CvnPtIx%2BYF%2FbMDMS1c4eBuLlwKrJJu08UFeaaFEhl0vYlnTAZUxQeqDUuKpgldJVSju%2FpA6Gq0ZjA%2BqlsN7iF3iZYkY63Ecc31P4QbEGCMGMzgnFWUxKFKBqDSPvqIjmIIB7tPJeBtilyKXMlB75X4jpHSFCRJpOSnPg1zsnxKrNMTRbOt5IG5s06kC1t47VG0X8R9uzNdQC%2Fygs8QmLnHHnbkoV6TR39QbIo3dhcpUXMD%2BJx1wj08OiRq3u31QEEEiAgYxjgUbFegMn3cFJDqkHA06%2BJSj%2Br%2FMvDTETO45GEqrYIxKqlAXpvh9mWkOK1nAjax%2BzlUHPGaVzhMHTMV8tZ0xvhJkV2BYSyRv%2BXLsKRMaj00k4RPW2eeOqvYaaLXt60hIjpCnDVu04qTt2KqzryQOIkw6%2BM8fRJQtCVVZlTbFMEwJZBo6peESVeQZEY9cjsMlAkqJJjmtDhSApzTSQDnTCeYYmS5aFADOfjgyiCh%2FIGkJtHW1GsrbLOoMqDn43sE4LBmtw1HH4pwZ8ohg%2F%2BBA1nC0gWCV%2B0rJENkKVaeGmuSE3VshYoY2gCxu9O9K7DD4x4PRH196IyKts1UFKSThypVXSR7uo6Me5n6I7AnrpvY8xiI6jtqlh0kVty7nPYQeRTg5LUn0XVS4yQU%2FqAMpfOFgPkFaZRu0IHiD7Y2zd31%2B698HX6RhcXMlI%2FjSWNyei%2FAYgqxIqPGhEubM7yOgZcl0SbfjlnAld0hX16YyM%2B%2BeSCK%2BUvtUltnqMKadgKlZj8Mi6whzvtAMEId2PuJv6gPqYXEbd4uTET3eSNRdKdqtczmjfKpdFsqQmq7n74RM7ulCNFgnNsln4g6UDiGtdsavyWtMackaYi7aHmdVFXAuTDpTNmSXHy07rx%2Bg6bOqZ1lUq2IeWCuzWszgGq9IFEVU0Sbvm0fLoEGdIVZdqLJU98aBap6n9f7OkloPqma4cpfWtZ6bfz%2F6JUFXP3CLuQo0i6setaRTnEewu%2B7UjyP4DiPZ62JUYxr2KdjHx9AYSJ4bx%2BI5axlLJH%2FIuAn6M6MyhkNpYteVKXS4lunKPS1Vd0pVV92wufAjZaq5NwbrlltNFHbXkgcY4BZk7NVtqXMCRbulT8ofQ8pwuVp1085wcy%2Btj8odeHqfnOCSPBj6bakQQyCNK0lEzdfUATKFJ8qQa3FWaCl9ooKnazlzRkQWMBgt0dzSw7AfubpD28MH9d5lGSyo26L4SyjJnDa8cOfnT3eSG68%2FplokqipSLYfgMzZLhxLzf0jsOt%2FNMlauzr%2F9h5I59wDQvopQHNC1pczZSzVHx%2FIBBmb%2FzVO1LsPADyprXKHhA9KkncVuYhlWFHYqdr5OtAQfzWAvDHNjCkBprI0bJMkiTYa6OB9crrYhHYh6u2Z62EIYNmEpxVutlcQ2iAKerCJ3Ts8y7d33%2BR7F4voW%2B5yN4hUJ2ETv%2F7ept3RSOIlnpum6s3WFGx3VjDX7dBHg20EVDejbOl4yLyRJwU7ZGMoWCgEsXUa53hsXiR2BnQ%2B2Ijg4RdRZvuOjVITrY7Xn7m6%2F5Hsrte%2FFOZb52zFlt34o3FJgcI2rNDFSwZfE%2Frv6znk5uF7AtUMgvnCRzP1yRyydzFBB9n0xue4%2BMfsbQ9yjjse4zgLPGN%2B4ySakHXBgpTQock84JdJuTyUwyXs3Dl5uil0OR8yeCMnK%2Be3J76PIAwl%2Fol3zbUbmLYDeagl0cY%2BxLrtsf5Ust%2BPHp6svU%2B%2F71ypeX362x%2FGrTymoTcLhZ3xDf6D2cYVtljtl6vW7INQ7bKfzR5hm%2FXVwseei08Ix%2FBvE9uT1GC8y2bcFlFBDZm8H%2FU8Vf4wV5bbKzDz7jghPEXw%2BgZ%2Bp%2BedUwBSbzKT3zYinZq%2FJ50ZkYelcTo%2FdMjP9lKPc9MSqqwiDbPAGdu5s%2FadxCkWbc48XcAqGLYwhCQgR2DqXu5KiyociV5Kiceqx3lhG6NHoxMqKelqEqRu8yQoBdXjx0UXuFLp0D6VafANV%2Bm%2BQDT7KQzHqvRshvn3Y5jJX8bCNk1wwYkKBb5BI95sQswltPhslPUb1OLyZ%2FfWGZmD6LALDKK9Ozz6MijnYKLR%2BXe7y7LXCMfJyj2QLVdHldbaZ5n1LLO%2BYSXuGbvyK0GH%2B5tVPbvTLGr1Y27CMKrM6i4ETQv1t0Ie%2F2XtGFo2%2BMu0UeIlNB7Q2Jht4X5dhjI8YwvG1x0%2B9nTJhC9MS8KbPKGOqhB%2BgkYTSh2QvsufchfgxHlcyFzOW%2Fy2hPvGWLxyToRnMqwwtVsm0T2147pga2Aarq3dF1gUS3T1rtxmrmWDYId9zdToV0eb27nTaGmF9iKZRC7FZDqmyBbQ4Av7aaHhlJqjmao2yfE%2FBDloVO2kK2v7VIqMs2P4l2RWUB3FJ4lwDpKJlslJivtgbGYDs26q84x4ni54Kq5sLwuXm01CmzT8Nlv%2FD5aCfD5cSpU3LuJnphxk3R715ypz6EXVKnwCMEKdiiOXsyYtim2NE6lfrdLz%2BKqqUO%2BjE4QrDI1CfRjg%2FUsGCVVbNU1fKGoRLRtyc2nbD%2FCVhEAZInt%2FKE0Xgy%2BzXsfip79PMw%2BqunhOB6OYuOOk2WjmXM5dt0B5MG2kxup3srK5mEAxJlDliUBBnD2DynnaJoKtuKnaE0sl3kwDBdtcA0Ts5OfAiz7yTdIka2BeUwdoI70Ssu7LftIYU3RO5YDsayjwUKteGAwt8inz7%2FZM32KLZ0IlDYFgxjqdy3PEe9Egfjg6Q4InNfMgRBJMo9zHyvmXNjdJLoV55xTru5c%2FDLLtN04MEv06oZjFrfSS7i76i1KtyGXnte0fZPaEax6loztFSLd1oVBt8GXxI8wi3r5e8c8abGzvnWnfIiZqjKxwOwBsFPtbxZpckYmoAxDuFEF3JBs2JY9bMGjbnore5Dg3aGoO6DyJt%2BiK%2FLCWn3nNiiUYesJedt8%2BNjjLAie4YYMNRedBL6V%2BmTHOvZ6R2a7gkST1yWt8FdstJxaWq7F%2Fjrxl6nDG9t5O%2FXG90SkuUFBre46%2F1ya4D9tcWwNhLj1QZxBtuxHqNLR9ExVkcdc4giRUJZ1mcpFj8LIFzimvdtkOlOQ3d0HYM9FVHWphACGUdiz%2FaKJ4Vz%2FkNK6OmWN%2Bev%2FfZsQpcJK2fDua3YD%2F5qTEwld0f0ZmM2XfsXOM32XwLq0cvS6GDukBYlh7BckFc0aargsyDiSTuAQ1mYHd2rP3nA2dFCYnX1Fb99FKRw8k19Ys55OF6QgVLxLaUYBxtrrAwLEExHmz8LIv4ECB3kpEgf7vAJEB1Ykqe5huMiCC1NtZCjqNBSZeQYEvSssaJ4qmYrnm5KrgFUyUGqqUPT0xTL8ZBlWq4DHQTlyktO9AmQraPl1mnLZfkc7%2FBdtRPlETxj8olSSDrOvSbpLvSA4RhEhzsKci3X8lTbMWzb9TQNjSVVg55tkkPbdFxNt02gOLqKINAN24KGAh1XhrbS29xvGGxp6ptX7TrzG%2FLr18n7QohYf0Ap%2F5T2LdPiCc%2F1z1JUy9n%2FWSl%2F9HPJk3KKMkV1afasb6a065zWiNgAEidrOSZaDlW2J5moOyu8LlEwIRUVgYVxLCfrJujwip2sYrJILwV%2FHhJHWl1xpNUzkOQlNFnI6RMt%2BlnJO3BxtGLDSO6ZvJxl4nRU2izAG3OXD%2B3WkRIQFn5CPckgRHhJj3nXWZizSERQhKkI7XP0MjIRNEWrVUUXFf7sfSNunjzVV%2B7STtYmOfgTxX72jbReqsSJaShIYd8kYk6fwr6x25307zIOVtMYuPd0gMI1U5reNX9kRymvE3Axtnek%2FQ7rzdZr641YYaLcH8E%2BSdU8SBEMQcWrAWUGDqra4qYiituz%2FrS%2B1pG40Emve4JeUqGTZ836MIot5sc9F1tkueLlRPH6B4kY7KnVWZEAfUqeZh7hxOem9BERYKkAS9HZTRmo7SzyMmBfvRKjlhfuG1KNtQFtYjxOvnpTYFUKM%2B4rvdpKbHWQXpVKgadXUXKfM%2F6yVNQzJjn7us7pPAG81Cb9rMy3q6oKoG2jSoGXZZLta2x%2BbuYU4p93xGVd3awAthSyO74CaEh7Ade0KoB6mU3N0odYZvPl%2B%2FgOU2Wvs9evv6V9fd51aTc%2BNHXilX19vntlvRfk0mssbVk67dKmMBrjtGw8kMHNP2OK6NX3%2Fwc%3D).

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
