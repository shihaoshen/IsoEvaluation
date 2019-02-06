#!/bin/sh
#
#___INFO__MARK_BEGIN__
##########################################################################
#
#  The Contents of this file are made available subject to the terms of
#  the Sun Industry Standards Source License Version 1.2
#
#  Sun Microsystems Inc., March, 2001
#
#
#  Sun Industry Standards Source License Version 1.2
#  =================================================
#  The contents of this file are subject to the Sun Industry Standards
#  Source License Version 1.2 (the "License"); You may not use this file
#  except in compliance with the License. You may obtain a copy of the
#  License at http://gridengine.sunsource.net/Gridengine_SISSL_license.html
#
#  Software provided under this License is provided on an "AS IS" basis,
#  WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING,
#  WITHOUT LIMITATION, WARRANTIES THAT THE SOFTWARE IS FREE OF DEFECTS,
#  MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE, OR NON-INFRINGING.
#  See the License for the specific provisions governing your rights and
#  obligations concerning the Software.
#
#  The Initial Developer of the Original Code is: Sun Microsystems, Inc.
#
#  Copyright: 2001 by Sun Microsystems, Inc.
#
#  All Rights Reserved.
#
##########################################################################
#___INFO__MARK_END__

# -- our name ---
#$ -N isoall
#$ -S /bin/sh

/bin/echo iso processing.

PATH=$PATH:$HOME/.local/bin:$HOME/bin
PATH=$PATH:/home/shens/bin/rMATS-ISO/lr2rmats/bin
export PATH

module load python/3.6.3
module load R/3.5.0
export PYTHONPATH=$HOME/py-ve/python3.6/site-packages
export R_LIBS=$HOME/R/R3.5/library

gtf=/mnt/isilon/xing_lab/shens/Annotation/Homo_sapiens.GRCh37.75.gtf
iso=/home/shens/bin/rMATS-ISO/
cd /mnt/isilon/xing_lab/shens/Run/Iso/CCLE
python /home/shens/bin/rMATS-ISO/rMATS-ISO.py --in-gtf $gtf --in-bam /mnt/isilon/xing_lab/shens/CCLE/bam/bam_input.ccle.list -o ./output_ccle/

qstat -ext
