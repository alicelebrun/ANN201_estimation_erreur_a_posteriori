function write_field(fname, name, field, Numtri, Coorneu)
  Nbpt = size(Coorneu,1);
  Nbtri = size(Numtri,1);

  fID = fopen(fname, 'w');
  fprintf(fID,'View "%s" {\n', name);

  for l=1:Nbtri
    S1=Coorneu(Numtri(l,1),:);
    S2=Coorneu(Numtri(l,2),:);
    S3=Coorneu(Numtri(l,3),:);

    fprintf(fID, ...
            '  ST (%s, %s, %s) {%s};\n', ...
            sprintf('%f,%f,0.0', S1), ...
            sprintf('%f,%f,0.0', S2), ...
            sprintf('%f,%f,0.0', S3), ...
            sprintf('%f, %f, %f', field(Numtri(l,:))));
  end
  fprintf(fID, '};\n');
  fclose(fID);
end
