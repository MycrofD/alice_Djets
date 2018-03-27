#!/bin/bash


root -l -b -q BkgMatrixRanCones.C'("fOutDzeroEvt.root","hDeltaPt_ptleadbin5_excluding","Djet5Excl",15)'
root -l -b -q BkgMatrixRanCones.C'("fOutDzeroEvt.root","hDeltaPt_ptleadbin10_excluding","Djet10Excl",15)'
root -l -b -q BkgMatrixRanCones.C'("fOutDzeroEvt.root","hDeltaPt_ptleadbin5_trans","Djet5Perp",15)'
root -l -b -q BkgMatrixRanCones.C'("fOutDzeroEvt.root","hDeltaPt_ptleadbin10_trans","Djet10Perp",15)'

root -l -b -q BkgMatrixRanCones.C'("fOutDzeroLeadEvt.root","hDeltaPt_ptleadbin5_excluding","Djet5Excl_Dlead",15)'
root -l -b -q BkgMatrixRanCones.C'("fOutDzeroLeadEvt.root","hDeltaPt_ptleadbin10_excluding","Djet10Excl_Dlead",15)'
root -l -b -q BkgMatrixRanCones.C'("fOutDzeroLeadEvt.root","hDeltaPt_ptleadbin5_trans","Djet5Perp_Dlead",15)'
root -l -b -q BkgMatrixRanCones.C'("fOutDzeroLeadEvt.root","hDeltaPt_ptleadbin10_trans","Djet10Perp_Dlead",15)'

root -l -b -q BkgMatrixRanCones.C'("fOutAllEvt.root","hDeltaPt_ptleadbin5_excluding","Alljet5Excl",25)'
root -l -b -q BkgMatrixRanCones.C'("fOutAllEvt.root","hDeltaPt_ptleadbin10_excluding","Alljet10Excl",25)'
root -l -b -q BkgMatrixRanCones.C'("fOutAllEvt.root","hDeltaPt_ptleadbin5_trans","Alljet5Perp",25)'
root -l -b -q BkgMatrixRanCones.C'("fOutAllEvt.root","hDeltaPt_ptleadbin10_trans","Alljet10Perp",25)'



