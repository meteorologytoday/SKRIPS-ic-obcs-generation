function o = read_json(filename)

    fid = fopen(filename);

    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);

    o = jsondecode(str);

end
