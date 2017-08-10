% Find the range of raw data
function result = chopnum(fin, sraw)
    if(nargin<2)
        sraw = 0;
    end
    hin = fopen(fin, 'r');
    rp = 0;
    result = [sraw sraw]; % start and end line to read
    status = 1;
    while(~feof(hin))
        tline = fgetl(hin);
        if(rp >= sraw)
            if(regexpi(tline, '^[^0-9+-.\s]'))
                if(2==status)
                    break;
                end
            else
                if(1==status)
                    result(1) = rp;
                    status = 2;
                end
            end
        end
        rp = rp + 1;
    end
    result(2) = rp - 1;
    fclose(hin);
end