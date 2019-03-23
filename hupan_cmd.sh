#!/bin/bash

IFS=' '
complete -W "qualSta trim alignRead sam2bam bamSta assemble alignContig extractSeq assemSta getUnalnCtg rmRedundant pTpG geneCov geneExist subSample gFamExist bam2bed fastaSta sim rmCtm blastAlign simSeq splitSeq getTaxClass genePre mergeNovGene" eupan

complete -W "qualSta mergeQualSta trim alignRead sam2bam bamSta assemble alignContig extractSeq assemSta mergeAssemSta getUnalnCtg mergeUnalnCtg rmRedundant pTpG geneCov mergeGeneCov geneExist subSample gFamExist bam2bed fastaSta sim rmCtm blastAlign simSeq splitSeq getTaxClass genePre mergeNovGene" eupanLSF

complete -W "qualSta mergeQualSta trim alignRead sam2bam bamSta assemble alignContig extractSeq assemSta mergeAssemSta getUnalnCtg mergeUnalnCtg rmRedundant pTpG geneCov mergeGeneCov geneExist subSample gFamExist bam2bed fastaSta sim rmCtm blastAlign simSeq splitSeq getTaxClass genePre mergeNovGene" eupanSLURM


