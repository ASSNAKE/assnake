include: results + 'metaphlan2/metaphlan2.py'
include: results + 'megahit/megahit_cross.py'
include: results + 'bwa/bwa.py'
include: results + 'centrifuge/centrifuge.py'
include: results + 'anvio/anvio.py'
include: results + 'trimmomatic/trimmomatic.py'
include: results + 'count/count.py'
include: results + 'fastqc/fastqc.py'

include: snakefiles + "bowtie2.py"
include: snakefiles + "megares.py"
include: snakefiles + "humann2.py"
include: snakefiles + "qiime2.py"

include: snakefiles + "ariba.py"
include: snakefiles + "prokka.py"

include: snakefiles + "general.py"
include: snakefiles + "preprocess.py"
include: snakefiles + "taxa.py"
include: snakefiles + "strain_finder.py"
include: snakefiles + "download.py"
include: snakefiles + "find_fungi.py"