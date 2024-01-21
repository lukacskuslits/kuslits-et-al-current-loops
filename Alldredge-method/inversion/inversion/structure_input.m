function param = structure_input(parname, loop_pars, ii, loop_nr, par_nr, nloops, fix)
%% Restructures the input loop parameters into the expected array format 
%% based on which ones of them are fixed
% Input parameters:
%------------------
% parname: name of the loop parameter which is to be estimated
% loop_pars: array containing the other loop parameters
% loop_nr: number of the loop (column of array 'param') for which the estimation is carried out (0 when
% it is performed on every loop simultaneously)
% par_nr: number corresponding to the position of the specified loop parameter in the input array
% nloops: number of current loops in the model
% fix: set containing the names of fixed loop parameters e.g.
% {'K', 'alpha', 'theta', 'phi'}
% Ouptut values:
%------------------
% param: restructured loop parameter array

%% HUN
% Az Alldredge-fele hurokparameterek hogyan erkeznek be attol fuggoen, melyeket rogzitjuk 
% Bemeno parameterek:
%--------------------
% parname: a becsleshez elokeszitendo hurokparameter neve
% loop_pars: tomb, mely a tobbi hurokparametert tartalmazza
% loop_nr: annak az aramhuroknak a szama (a 'param' tomb oszlopa), amelyre eppen vegezzuk a becslest (0 when
% ha minden hurokra egyszerre vegezzuk)
% par_nr: a megadott hurokparameter pozicioja a bemeno (hurok)parameterek
% tombjeben
% nloops: aramhurkok szama a modellben
% fix: azon hurokparameterek nevei, amelyeket rogzitunk a becsles soran, pl.
% {'K', 'alpha', 'theta', 'phi'}
% Eredmeny:
%--------------------
% param: atalakitott hurokparameter tomb
%%
if strcmp(fix, 'none')
    param=loop_pars(ii, par_nr);
elseif strcmp(fix, [parname,'_loop'])
       param0 = load([parname,'0.mat']);
       param0 = param0.([parname,'0']);
    if loop_nr == 1
           param = [loop_pars(ii), param0(loop_nr+1:nloop)];
    elseif loop_nr == nloops
           param=[param0(1:loop_nr), loop_pars(ii)];
    else
           param=[param0(1:loop_nr), loop_pars(ii), param0(loop_nr+1:nloop)];
    end
else
     param=loop_pars(ii);
end