function stn = get_IGOSS_mld(stn)
%function stn = get_IGOSS_mld(stn)
%
% Download Mixed-Layer Depth climatology from IGOSS

  rawdat = loaddap('http://iridl.ldeo.columbia.edu/SOURCES/.IGOSS/.sio/.climatology/.mld/Y/%2828N%29%2824N%29RANGEEDGES/X/%2883W%29%2877W%29RANGEEDGES/dods');

return;
