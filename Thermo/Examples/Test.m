clear,clf, addpath ../ ../EOS ../Solutions/ ../Utilities/
run_name = 'EF21_red';
T        = linspace(400,4000,51);
P        = linspace( 5 ,  35,51)*1e9;  % Pa
solmod   = 'solution_models_H18';
ngrid    = 4; % number of P-T grid refinements, each time, the P-T grid resolution is doubled.
eps_solv = 1;
Cname   = {'Mg' 'Si' 'O' };
Nsys     = [2    1    4];
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'fo,tc-ds633','mwd,tc-ds633','mrw,tc-ds633','mpv,tc-ds633','per,tc-ds633','Melt(H18)'};
td       = init_thermo(phs_name,Cname,solmod);
for i = 1:length(phs_name),td(i).dz(:) = 1/4;end
p         = props_generate(td);     % generate endmember proportions
refine_id = ones(length(T)*length(P),1);
% Minimization refinement
for i_grid = 1:ngrid
    [T2d,P2d] = ndgrid(T,P);
    parfor iPT = 1:length(T2d(:))
        if refine_id(iPT) == 1
            [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td);
            disp(iPT/length(T2d(:)))
        end
    end
end
clear A
for i = 1:length(pc_id_ref)
    A(i) = max(pc_id_ref{i});
end

pcolor(T,P/1e9,reshape(A,length(T),[])');
shading interp
save_the_image('tiff',300,[14,14],'Fig')


