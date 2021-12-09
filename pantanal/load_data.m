RivMAP_path = 'C:\Users\renan\Unifesp\Iniciação Científica\ESTIMATIVA DA CONCENTRAÇÃO DE SEDIMENTOS EM SUSPENSÃO NA BACIA HIDROGRÁFICA DO RIO TAQUARIMS\rivmap_pantanal\pantanal';
cd(RivMAP_path)

from = 2013;
to = 2021;


for year = from:to
    index = year-from + 1;

    img_name= strcat('rivmapcuiaba_menor/cuiaba_menor_active_channel_binary_mask_', string(year), '.tif');
    img_number = imread(img_name);
    img = logical(img_number);
    
    meta = struct();
    meta.year = year;
    meta.exit_sides = 'EW';
    meta.Wn = 10;
    rivmapcuiaba_menor(index).meta = meta;
    
    im = struct();
    im.st = img;
    im.hc = img;
    rivmapcuiaba_menor(index).im = im;
     
end

save('rivmapcuiaba_menor');