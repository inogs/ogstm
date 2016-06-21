#! /bin/bash

module load autoload nco

TEST='TEST01'
MODDIR="../${TEST}/wrkdir/MODEL"

pw=$PWD

cd ${MODDIR}/AVE_FREQ_1/

# looping on timesteps to merge all variables
  for var in ABIO_eps Acae Acan ALK CO2 CO3 cxoO2 D1m D2m D6m D7m D8m D9m DIA_B_1 DICae DICan EICE eiPPY_iiP1 eiPPY_iiP2 eiPPY_iiP3 eiPPY_iiP4 ELiPPY_iiP1 ELiPPY_iiP2 ELiPPY_iiP3 ELiPPY_iiP4 eO2mO2 EPR exPPYc exPPYs flPTN6r fP1Z4c fR2B1c G13c G13h G23c G23h G2o G3h HCO3 jbotB1c jbotB1n jbotB1p jbotO3h jbotP1c jbotP1n jbotP1p jbotP1l jbotP1s  jbotP2c jbotP2l jbotP2n jbotP2p jbotP3c jbotP3l jbotP3n jbotP3p jbotP4c jbotP4l jbotP4n jbotP4p jbotZ3c jbotZ3n jbotZ3p jbotZ4c jbotZ4n jbotZ4p jbotZ5c jbotZ5n jbotZ5p jbotZ6c jbotZ6n jbotZ6p jG2K3o jG2K7o jK13K3n jK15K5s jK31K21p jK34K24n jK36K26r jP1Y3c jP2Y3c jP3Y3c jP4Y3c jRIY3c jRIY3n jRIY3p jRIY3s jsurB1c jsurB1n jsurB1p jsurO3h jsurP1c jsurP1n jsurP1p jsurP1l jsurP1s  jsurP2c jsurP2l jsurP2n jsurP2p jsurP3c jsurP3l jsurP3n jsurP3p jsurP4c jsurP4l jsurP4n jsurP4p jsurZ3c jsurZ3n jsurZ3p jsurZ4c jsurZ4n jsurZ4p jsurZ5c jsurZ5n jsurZ5p jsurZ6c jsurZ6n jsurZ6p jZIY3c K11p K14n K16r K1p K21p K24n K26r K3n K4n K5s K6r M11p M14n M1p M21p M24n M3n M4n M5s M6r O3h_Ben O3h OArag OCalc pCO2ae pCO2an pCO2 pHae pHan pH ppP1_Benc ppP1_Benn ppP1_Benp ppP1_Bens ppP2_Benc ppP2_Benn ppP2_Benp ppP2_Bens ppP3_Benc ppP3_Benn ppP3_Benp ppP3_Bens ppP4_Benc ppP4_Benn ppP4_Benp ppP4_Bens ppreH1 ppreH2 ppruH1 ppruH2 Q11c Q11n Q11p qlcPPY_iiP1 qlcPPY_iiP2 qlcPPY_iiP3 qlcPPY_iiP4 qlcPPY_iiZ3 qlcPPY_iiZ4 qlcPPY_iiZ5 qlcPPY_iiZ6 qncOMT_iiR1 qncOMT_iiR2 qncOMT_iiR3 qncOMT_iiR6 qncPBA_iiB1 qncPPY_iiP1 qncPPY_iiP2 qncPPY_iiP3 qncPPY_iiP4 qpcMEZ_iiZ3 qpcMEZ_iiZ4 qpcMIZ_iiZ5 qpcMIZ_iiZ6 qpcOMT_iiR1 qpcOMT_iiR2 qpcOMT_iiR3 qpcOMT_iiR6 qpcPBA_iiB1 qpcPPY_iiP1 qpcPPY_iiP2 qpcPPY_iiP3 qpcPPY_iiP4 qscOMT_iiR1 qscOMT_iiR2 qscOMT_iiR3 qscOMT_iiR6 qscPPY_iiP1 qscPPY_iiP2 qscPPY_iiP3 qscPPY_iiP4 reATn reATp reBn reBTn reBTp rePTn rePTp resPBAc  resPPYc resZOOc RI_Fc RI_Fn RI_Fp RI_Fs rrATo ruBn ruPPYc ruPPYn ruPPYp ruPPYs ruZOOc sediMEZ_iiZ3 sediMEZ_iiZ4 sediMEZ_iiZ5 sediMEZ_iiZ6 shiftD1m shiftD2m sunPPY_iiP1 sunPPY_iiP2 sunPPY_iiP3 sunPPY_iiP4 SUNQ ThereIsLight xEPS Y1c Y1n Y1p Y2c Y2n Y2p Y3c Y3n Y3p Y4c Y4n Y4p Y5c Y5n Y5p ZI_Fc ZI_Fn
  do
  rm ./*.$var'.nc'
  done

for timesteps in *.O2o.nc; 
  do
  timestep=$(echo $timesteps| cut -b -21)
  echo $timestep
  for var in $timestep*.nc;
    do
    ncks -Ah  $var 'F1_'$timestep.nc
    done
  done

cd ../AVE_FREQ_2/
  for var in ABIO_eps Acae Acan ALK CO2 CO3 cxoO2 D1m D2m D6m D7m D8m D9m DIA_B_1 DICae DICan EICE eiPPY_iiP1 eiPPY_iiP2 eiPPY_iiP3 eiPPY_iiP4 ELiPPY_iiP1 ELiPPY_iiP2 ELiPPY_iiP3 ELiPPY_iiP4 eO2mO2 EPR exPPYc exPPYs flPTN6r fP1Z4c fR2B1c G13c G13h G23c G23h G2o G3h HCO3 jbotB1c jbotB1n jbotB1p jbotO3h jbotP1c jbotP1n jbotP1p jbotP1l jbotP1s  jbotP2c jbotP2l jbotP2n jbotP2p jbotP3c jbotP3l jbotP3n jbotP3p jbotP4c jbotP4l jbotP4n jbotP4p jbotZ3c jbotZ3n jbotZ3p jbotZ4c jbotZ4n jbotZ4p jbotZ5c jbotZ5n jbotZ5p jbotZ6c jbotZ6n jbotZ6p jG2K3o jG2K7o jK13K3n jK15K5s jK31K21p jK34K24n jK36K26r jP1Y3c jP2Y3c jP3Y3c jP4Y3c jRIY3c jRIY3n jRIY3p jRIY3s jsurB1c jsurB1n jsurB1p jsurO3h jsurP1c jsurP1n jsurP1p jsurP1l jsurP1s  jsurP2c jsurP2l jsurP2n jsurP2p jsurP3c jsurP3l jsurP3n jsurP3p jsurP4c jsurP4l jsurP4n jsurP4p jsurZ3c jsurZ3n jsurZ3p jsurZ4c jsurZ4n jsurZ4p jsurZ5c jsurZ5n jsurZ5p jsurZ6c jsurZ6n jsurZ6p jZIY3c K11p K14n K16r K1p K21p K24n K26r K3n K4n K5s K6r M11p M14n M1p M21p M24n M3n M4n M5s M6r O3h_Ben O3h OArag OCalc pCO2ae pCO2an pCO2 pHae pHan pH ppP1_Benc ppP1_Benn ppP1_Benp ppP1_Bens ppP2_Benc ppP2_Benn ppP2_Benp ppP2_Bens ppP3_Benc ppP3_Benn ppP3_Benp ppP3_Bens ppP4_Benc ppP4_Benn ppP4_Benp ppP4_Bens ppreH1 ppreH2 ppruH1 ppruH2 Q11c Q11n Q11p qlcPPY_iiP1 qlcPPY_iiP2 qlcPPY_iiP3 qlcPPY_iiP4 qlcPPY_iiZ3 qlcPPY_iiZ4 qlcPPY_iiZ5 qlcPPY_iiZ6 qncOMT_iiR1 qncOMT_iiR2 qncOMT_iiR3 qncOMT_iiR6 qncPBA_iiB1 qncPPY_iiP1 qncPPY_iiP2 qncPPY_iiP3 qncPPY_iiP4 qpcMEZ_iiZ3 qpcMEZ_iiZ4 qpcMIZ_iiZ5 qpcMIZ_iiZ6 qpcOMT_iiR1 qpcOMT_iiR2 qpcOMT_iiR3 qpcOMT_iiR6 qpcPBA_iiB1 qpcPPY_iiP1 qpcPPY_iiP2 qpcPPY_iiP3 qpcPPY_iiP4 qscOMT_iiR1 qscOMT_iiR2 qscOMT_iiR3 qscOMT_iiR6 qscPPY_iiP1 qscPPY_iiP2 qscPPY_iiP3 qscPPY_iiP4 reATn reATp reBn reBTn reBTp rePTn rePTp resPBAc  resPPYc resZOOc RI_Fc RI_Fn RI_Fp RI_Fs rrATo ruBn ruPPYc ruPPYn ruPPYp ruPPYs ruZOOc sediMEZ_iiZ3 sediMEZ_iiZ4 sediMEZ_iiZ5 sediMEZ_iiZ6 shiftD1m shiftD2m sunPPY_iiP1 sunPPY_iiP2 sunPPY_iiP3 sunPPY_iiP4 SUNQ ThereIsLight xEPS Y1c Y1n Y1p Y2c Y2n Y2p Y3c Y3n Y3p Y4c Y4n Y4p Y5c Y5n Y5p ZI_Fc ZI_Fn
  do
  rm ./*.$var'.nc'
  done

# looping on timesteps to merge all variables                                                                                                                                                                    
 for timesteps in *.O2o.nc;
  do
  timestep=$(echo $timesteps| cut -b -21)
  echo $timestep
  for var in $timestep*.nc;
    do
    ncks -Ah  $var 'F2_'$timestep.nc
    done
  done

cd $pw

mv ${MODDIR}/AVE_FREQ_1/F1_* .
mv ${MODDIR}/AVE_FREQ_2/F2_* .


#ncrcat -O -h ${MODDIR}/AVE_FREQ_1/ave*${VAR}.nc aux.nc
#ncks -O --mk_rec_dmn time aux.nc ave.${VAR}.nc
#rm aux.nc


