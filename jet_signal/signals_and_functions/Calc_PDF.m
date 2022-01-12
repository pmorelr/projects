function [Bin,PDF]=Calc_PDF(SIG,NbBin)
% fonction qui calcule la densit� de probabilit� (PDF) d'un signal SIG
% donn� sur un nombre d'intervals NbBin
%
% calcul de l'histogramme:
[PDF,Bin]=hist(SIG,NbBin); 
% normalisation pour avoir la densit� de probabilit�:
Int=cumtrapz(Bin,PDF);
Int=Int(end);
PDF=PDF/Int;
