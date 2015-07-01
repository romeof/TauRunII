#!/bin/bash
#Specify needed variables
#int
#varType=int
#varTYPE=INT
#capLetter=I
#varList=(rtau_pt rtau_eta rtau_phi rtau_en)
#varLast=jet_chtrksnpvtt
#varNum=jn
#varCount=jn
#double
varType=double
varTYPE=DOUBLE
capLetter=D
varList=(pv_avfbs_nchi2 unbpv_avfbs_nchi2 diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2)
#pvsv_dist3d_val pvsv_dist2d_val pvsv_dist3d_sig pvsv_dist2d_sig)
#recojettau_pt recojettau_eta recojettau_phi recojettau_en)
#unbpv_KVF_nv unbpv_KVFbs_nv unbpv_AVF_nv unbpv_AVFbs_nv)
#vtxKVF_x vtxKVF_y vtxKVF_z vtxKVFbs_x vtxKVFbs_y vtxKVFbs_z vtxAVF_x vtxAVF_y vtxAVF_z vtxAVFbs_x vtxAVFbs_y vtxAVFbs_z)
#genditauvtx_x genditauvtx_y genditauvtx_z)
#pvtauvtx_genditauvtx_x pvtauvtx_genditauvtx_y pvtauvtx_genditauvtx_z unbpv_genditauvtx_x unbpv_genditauvtx_y unbpv_genditauvtx_z)
varLast=
varNum=ele_num
varCount=gl
#Print info
echo " "
#Declare variables
echo -e " $varType \c"
pos=0
for count in ${varList[@]}; 
do
  if [ "${varList[$pos]}" != "$varLast" ] 
  then
   echo -e "${varList[$pos]}, \c"
  else
   echo "${varList[$pos]};"
  fi
  let pos=pos+1
done
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  echo "  ${varList[$pos]} = DEF_VAL_$varTYPE;"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->Branch(\"${varList[$pos]}\", &${varList[$pos]}, \"${varList[$pos]}/$capLetter\");"
  let pos=pos+1
done
echo " "
#Set branch address
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->SetBranchAddress(\"${varList[$pos]}\", &${varList[$pos]});"
  let pos=pos+1
done
echo " "
#Analyzer
pos=0
for count in ${varList[@]}; 
do
  echo "tree->${varList[$pos]} = ;"
  let pos=pos+1
done
echo " "
