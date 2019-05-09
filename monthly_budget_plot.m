1;

RAPFX = 'erai';
TURPFX = ['ndbc_' RAPFX '_30a'];
q0fld = [TURPFX '_net_heat_flux'];
qtfld = [q0fld '_term'];
QEPFX = 'ww3_fkeys_qe';
dTfld = [TURPFX '_' QEPFX '_dt'];


absfld = 'ndbc_sea_t';
accfld = qtfld;

nboots=30;

if ( ~exist('stats','var') || ~isfield(stats,absfld) )
  disp(absfld);
  [stats.(absfld).m,stats.(absfld).s] = bootmon(stn,absfld,nboots);
end;

if ( ~isfield(stats,accfld) )
  disp(accfld);
  [stats.(accfld).m,stats.(accfld).s] = bootmon(stn,accfld,nboots);
end;

plotbootmon(stats,absfld,accfld);
