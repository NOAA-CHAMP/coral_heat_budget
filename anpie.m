1;

[qlh,qsh,qlw,gqsw,qb,adv,dif,hc] = ...
    intersect_tses(stn.ndbc_erai_erai_30a_latent_flux,...
    stn.ndbc_erai_erai_30a_sensible_flux,...
    stn.erai_ndbc_lrf,...
    stn.absorbed_erai_ndbc_srf,...
    stn.b_erai_qbo,...
    stn.fq_erai_avhrr_advected_heat_flux,...
    stn.avhrr_weekly_diffused_heat_flux,...
    stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc_flux);

alldat = [qlh.data,qsh.data,qlw.data,gqsw.data,qb.data,adv.data,dif.data,hc.data];

%cf=@iqr;
%cf=@nanstd;
%cf=@(x)(abs(nanmedian(x)));
%cf=@(x)(abs(nanmean(x)));
%cf=@(x)(nanmedian(abs(x)));
cf=@(x)(nanmean(abs(x)));

dat = [cf(qlh.data),cf(qsh.data),cf(qlw.data),cf(gqsw.data),cf(qb.data),cf(adv.data),cf(dif.data),cf(hc.data)];

legs = {'Q_L_H','Q_S_H','Q_L_W','\gammaQ_S_W','Q_b','F_qu\bullet\nabla_hT_s','K_H\nabla_h^2T_s','u_H_C\bullet\nablaT_s(\beta,Q,h)'};

fmg;
pie(dat,[0,0,0,0,1,0,1,0]);
legend(legs,...
       'Location','Best');
axis square
colormap(gray);
