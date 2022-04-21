
clc
clear all
close all

load('essaisLaura')
% force en kN et déplacement en mm

essais = {};

figure
hold on

dec=0;
for e=1:20
    
    %% Récupéraiton des données

    colonne = e*3-2;
    
    % récupère le temps
    coTemps = colonne+dec;
    temps = table2array(essaisLaura(:,coTemps));

    % récupère le déplacement
    if e == 9
        % l'essai 9 à 2 colonnes de déplacement on prend la 2 eme colonne
        % Que c'est til passé ?
        dec=1;
    end
    coDeplacement = colonne+1+dec;
    deplacement = table2array(essaisLaura(:,coDeplacement));

    % récupère la force
    coForce = colonne+2+dec;
    force = table2array(essaisLaura(:,coForce));

    %enlève les Nan
    temps=temps(~isnan(temps));
    deplacement=deplacement(~isnan(deplacement));
    force=force(~isnan(force));
    
    % mets les vecteurs à la bonne dimension
    dim = min([length(temps),length(deplacement),length(force)]);
    liste = 1:dim;
    
    % stocke les résultats
    temps = temps(liste);
    deplacement = deplacement(liste);
    force = force(liste);
    force = smoothdata(force,'movmean',10); % 'movmean','movmedian','gaussian','lowess','loess','rlowess','rloess','sgolay' 
    
    result = [temps, deplacement, force];
    essais{e} = result;

    %% On determine la force critique
    
    % detecte a partir de quand le force ne décroit plus
    % Pour ça il faut detecter quand on a le dernier 0
    debut = min(find(force>0.05));

    incTemps = temps(debut+1:end)-temps(debut:end-1);
    derive = (force(debut+1:end)-force(debut:end-1))./diff(deplacement(debut:end));
    
    incMax = min(find(derive<1e-12))+debut;

    forceMax(e) = force(incMax);
    depMax(e) = deplacement(incMax);

    % figure force déplacement
    plot(deplacement, force)    
    legends(e) = "Essai "+e+"";

end

xlim([0,10])
ylim([0,4.5])
grid on
xlabel('Déplacement en mm')
ylabel('Force de compression en kN')
plot(depMax, forceMax, 'rx', LineWidth=1.5)
legend(legends)

fprintf('fc [%3.2f; %3.2f]\n', min(forceMax), max(forceMax))

fMean = mean(forceMax)
fEcartType = std(forceMax)
