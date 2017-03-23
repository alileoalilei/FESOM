#!/bin/bash
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!! SUBMIT ALL POSTPROC ROUTINES AT ONCE !!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

date

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!! SUBMIT TIMESERIES DIAGNOSTICS  !!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#bsub -K < SBMT_FEOM_DIAG.sh
#bsub -K < SBMT_FEOM_MARM.sh


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!! SUBMIT THALWEG DIAGNOSTICS !!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#TMPLFILE=SBMT_FEOM_THAL_TeMPLaTe.sh
#SBMTFILE=SBMT_FEOM_THAL.sh
#
#INITIALDAY=360; FINALDAY=360; EXPDEF=BAS; EXPNUM=01; YEAR=2008
#sed -e 's;INIDAY;'${INITIALDAY}';g' -e \
#       's;ENDDAY;'${FINALDAY}';g' -e \
#       's;EXPNAME;'${EXPDEF}';g' -e \
#       's;EXPNO;'${EXPNUM}';g' -e \
#       's;EXPYEAR;'${YEAR}';g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}

#TMPLFILE=SBMT_FEOM_THAL_MEAN_TeMPLaTe.sh
#SBMTFILE=SBMT_FEOM_THAL_MEAN.sh
#
#FSTMON=293; LSTMON=293; EXPDEF=BLK; EXPNUM=01; YEAR=2011
#sed -e 's;INIMON;'${FSTMON}';g' -e \
#       's;ENDMON;'${LSTMON}';g' -e \
#       's;EXPNAME;'${EXPDEF}';g' -e \
#       's;EXPNO;'${EXPNUM}';g' -e \
#       's;EXPYEAR;'${YEAR}';g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}


#TMPLFILE=SBMT_FEOM_THAL_ANNUAL_MEAN_TeMPLaTe.sh
#SBMTFILE=SBMT_FEOM_THAL_MEAN.sh
#
#FSTMON=1; LSTMON=1; EXPDEF=BLK; EXPNUM=01; YEAR=2008
#sed -e 's;INIMON;'${FSTMON}';g' -e \
#       's;ENDMON;'${LSTMON}';g' -e \
#       's;EXPNAME;'${EXPDEF}';g' -e \
#       's;EXPNO;'${EXPNUM}';g' -e \
#       's;EXPYEAR;'${YEAR}';g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!! SUBMIT SECTION DIAGNOSTICS !!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#TMPLFILE=SBMT_FEOM_SECT_TeMPLaTe.sh
#SBMTFILE=SBMT_FEOM_SECT.sh
#
#INITIALDAY=1; FINALDAY=1; EXPDEF=BLK; EXPNUM=02; YEAR=2008
#sed -e 's;INIDAY;'${INITIALDAY}';g' -e \
#       's;ENDDAY;'${FINALDAY}';g' -e \
#       's;EXPNAME;'${EXPDEF}';g' -e \
#       's;EXPNO;'${EXPNUM}';g' -e \
#       's;EXPYEAR;'${YEAR}';g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}


#TMPLFILE=SBMT_FEOM_SECT_MEAN_TeMPLaTe.sh
#SBMTFILE=SBMT_FEOM_SECT_MEAN.sh
#
#FSTMON=1; LSTMON=1; EXPDEF=BLK; EXPNUM=02; YEAR=2008
#sed -e 's;INIMON;'${FSTMON}';g' -e \
#       's;ENDMON;'${LSTMON}';g' -e \
#       's;EXPNAME;'${EXPDEF}';g' -e \
#       's;EXPNO;'${EXPNUM}';g' -e \
#       's;EXPYEAR;'${YEAR}';g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}
#
#
#TMPLFILE=SBMT_FEOM_SECT_ANNUAL_MEAN_TeMPLaTe.sh
#SBMTFILE=SBMT_FEOM_SECT_MEAN.sh
#
#FSTMON=1; LSTMON=365; EXPDEF=BLK; EXPNUM=02; YEAR=2008
#sed -e 's;INIMON;'${FSTMON}';g' -e \
#       's;ENDMON;'${LSTMON}';g' -e \
#       's;EXPNAME;'${EXPDEF}';g' -e \
#       's;EXPNO;'${EXPNUM}';g' -e \
#       's;EXPYEAR;'${YEAR}';g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}
#
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!! SUBMIT FORCING DIAGNOSTICS !!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
TMPLFILE=SBMT_FEOM_FORC_ANNUAL_MEAN_TeMPLaTe.sh
SBMTFILE=SBMT_FEOM_FORC_MEAN.sh
#
## Variables are WSTC, WWRK, FORC and BUOY
#
FSTDAY=1; LSTDAY=365; EXPDEF=BLK; EXPNUM=02; YEAR=( 2008 2009 2010 )
sed -e 's;^INITIALDAY=.*$;INITIALDAY='${FSTDAY}';g' -e\
       's;^FINALDAY=.*$;FINALDAY='${LSTDAY}';g' -e\
       's;^EXPDEF=.*$;EXPDEF='${EXPDEF}';g' -e\
       's;^EXPNUM=.*$;EXPNUM='${EXPNUM}';g' -e\
       's;^VAR=.*$;VAR=(WWRK);g' -e\
       's;^YEARS=.*$;YEARS=( 2008 2009 2010 2011 2012 2013 );g' \
       ${TMPLFILE} > ${SBMTFILE}

bsub < ${SBMTFILE}
#
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!! SUBMIT OCEAN DIAGNOSTICS !!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#TMPLFILE=SBMT_FEOM_VORT_ANNUAL_MEAN_TeMPLaTe.sh
#SBMTFILE=SBMT_FEOM_VORT_MEAN.sh
#
#FSTDAY=1; LSTDAY=365; EXPDEF=BLK; EXPNUM=02; YEAR=( 2008 2009 2010 )
#sed -e 's;^INITIALDAY=.*$;INITIALDAY='${FSTDAY}';g' -e\
#       's;^FINALDAY=.*$;FINALDAY='${LSTDAY}';g' -e\
#       's;^EXPDEF=.*$;EXPDEF='${EXPDEF}';g' -e\
#       's;^EXPNUM=.*$;EXPNUM='${EXPNUM}';g' -e\
#       's;^YEARS=.*$;YEARS=(2008 2009 2010 2011 2012 2013);g' \
#       ${TMPLFILE} > ${SBMTFILE}
#
#bsub < ${SBMTFILE}
#
#TMPLFILE=SBMT_FEOM_FLUX_TeMPLaTe.sh
#SBMTFILE=SBMT_FEOM_FLUX.sh
#
#FSTDAY=1; LSTDAY=365; EXPDEF=BBI; EXPNUM=nc; YEAR=( 2008 2009 2010 )
#sed -e 's;^INITIALDAY=.*$;INITIALDAY='${FSTDAY}';g' -e\
#       's;^FINALDAY=.*$;FINALDAY='${LSTDAY}';g' -e\
#       's;^EXPDEF=.*$;EXPDEF='${EXPDEF}';g' -e\
#       's;^EXPNUM=.*$;EXPNUM='${EXPNUM}';g' -e\
#       's;^YEARS=.*$;YEARS=( 2008 );g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}
#

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!! SUBMIT ENSEMBLE DIAGNOSTICS !!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#TMPLFILE=SBMT_ENS_SECT_TeMPLaTe.sh
#SBMTFILE=SBMT_ENS_SECT.sh
#
#INITIALDAY=1 FINALDAY=4; EXPDEF=INI; EXPNUM=01; YEAR=2009
#sed -e 's;INIDAY;'${INITIALDAY}';g' -e \
#       's;ENDDAY;'${FINALDAY}';g' -e \
#       's;EXPNAME;'${EXPDEF}';g' -e \
#       's;EXPNO;'${EXPNUM}';g' -e \
#       's;EXPYEAR;'${YEAR}';g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!! SUBMIT OBSERVATION TOOLS   !!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#TMPLFILE=SBMT_CREATE_OBS_TeMPLaTe.sh
#SBMTFILE=SBMT_CREATE_OBS.sh
#
#INITIALDAY=1; FINALDAY=1; EXPDEF=WFX; EXPNUM=04; YEAR=2009
#sed -e 's;INIDAY;'${INITIALDAY}';g' -e \
#       's;ENDDAY;'${FINALDAY}';g' -e \
#       's;EXPNAME;'${EXPDEF}';g' -e \
#       's;EXPNO;'${EXPNUM}';g' -e \
#       's;EXPYEAR;'${YEAR}';g' \
#       ${TMPLFILE} > ${SBMTFILE}
#bsub < ${SBMTFILE}
date
