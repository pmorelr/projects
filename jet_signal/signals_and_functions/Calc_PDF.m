function [Bin,PDF]=Calc_PDF(SIG,NbBin)
% fonction qui calcule la densité de probabilité (PDF) d'un signal SIG
% donné sur un nombre d'intervals NbBin
%
% calcul de l'histogramme:
[PDF,Bin]=hist(SIG,NbBin); 
% normalisation pour avoir la densité de probabilité:
Int=cumtrapz(Bin,PDF);
Int=Int(end);
PDF=PDF/Int;
