function locStr=CWLocationString(pt_name,tar_loc),

locStr='';
%# Patients with non-standard DVH names
anom_pts = {'Cohen','Lapid','Menaker','Pisano','Sutton','Underkoffler'};

matches = strcmp(anom_pts,pt_name);

if sum(matches)>0,
    if tar_loc(2)=='U',
        locStr='CW2SUP';
    else
        locStr='CW2INF';
    end
else %# not anomoly, just determine L/R
     if tar_loc(1)=='R',
        locStr='CW2RT';
    else
        locStr='CW2LT';
     end
end
    
end
