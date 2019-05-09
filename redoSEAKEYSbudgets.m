1;

yOrN = questdlg('Are you SURE you want to rerun ALL? This may take anywhere from 30 minutes up to TWO DAYS, depending on options...');

if ( strcmpi(yOrN(1),'y') )

  closeAll = false;

  diary redoallbudgets.log
  disp(['STARTING ',upper(mfilename)]);
  timenow;

  % Meteorology reanalysis
  %for re={'ncep','erai'};
  for re={'erai'};
    % Km-scale (advection-diffusion)
    %for km={'gom_hycom','avhrr_weekly'};
    for km={'avhrr_weekly'};
      % Primary meteorology (wind, air temp)
      %for mt={'erai','ndbc'};
      for mt={'ndbc'};
        % Tides
        %for td={'tmd_tide','tpxo_tide'};
        for td={'tpxo_tide'};
          % Waves
          %for wv={'ndbc','ww3','erai'};
          for wv={'erai'};

            %for cst={'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1'};
            for cst={'fwyf1','mlrf1','lonf1','sanf1'};
              try, stn=optimize_station_heat_budget(cst{:},re{:},km{:},mt{:},td{:},wv{:});
              catch, catchwarn(cst{:}); end;
              stn=[]; clear stn; if (closeAll);close all;end; pause(1);
            end;
            pause(1);

            % Sombrero AVHRR SST gradients were just too bad to save
            try, stn=optimize_station_heat_budget('smkf1',re{:},'none',mt{:},td{:},wv{:});
            catch, catchwarn('smkf1'); end;
            stn=[]; clear stn; if (closeAll);close all;end; pause(1);
            % % Try merging CT and NDBC sea temps. for more complete record
            % try, stn=optimize_station_heat_budget('smkf1',re{:},km{:},mt{:},td{:},wv{:},'merged_sea_t');
            % catch, catchwarn('smkf1'); end;
            % stn=[]; clear stn; if (closeAll);close all;end; pause(1);

            % Dry Tortugas tidal currents from TPXO miniscule compared to GoM ("TMD") model
            try, stn=optimize_station_heat_budget('dryf1',re{:},km{:},mt{:},'tmd_tide',wv{:});
            catch, catchwarn('dryf1'); end;
            stn=[]; clear stn; if (closeAll);close all;end; pause(1);

            % Two different sea temperature records to process at Looe Key
            try, stn=optimize_station_heat_budget('looe1',re{:},km{:},mt{:},td{:},wv{:},'mc_seatemp');
            catch, warning('ERROR for looe1 MC'); e=lasterror; disp(e.identifier); disp(e.message); disp(e.stack(1)); end;
            stn=[]; clear stn; if (closeAll);close all;end; pause(1);
            try, stn=optimize_station_heat_budget('looe1',re{:},km{:},mt{:},td{:},wv{:},'ad_seatemp');
            catch, warning('ERROR for looe1 AD'); e=lasterror; disp(e.identifier); disp(e.message); disp(e.stack(1)); end;
            stn=[]; clear stn; if (closeAll);close all;end; pause(1);

          end;
        end;
      end;
    end;
  end;

  timenow;
  diary off

end;

clear cst km mt re td wv e yOrN closeAll
