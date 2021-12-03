RivMAP_path = 'C:\Users\renan\Unifesp\Iniciação Científica\ESTIMATIVA DA CONCENTRAÇÃO DE SEDIMENTOS EM SUSPENSÃO NA BACIA HIDROGRÁFICA DO RIO TAQUARIMS\rivmap_pantanal\pantanal';
cd(RivMAP_path)

from = 2013;
to = 2021;


for year = from:to
    index = year-from + 1;

    img_name= strcat('rivmapsao_lourenco/sao_lourenco_active_channel_binary_mask_', string(year), '.tif');
    img_number = imread(img_name);
    img = logical(img_number);
    
    rivmapsao_lourenco(index).meta.year = year
    rivmapsao_lourenco(index).meta.exit_sides = 'NS';
    rivmapsao_lourenco(index).meta.Wn = 30;
    
    rivmapsao_lourenco(index).im = struct();
    rivmapsao_lourenco(index).im.st = img;
    rivmapsao_lourenco(index).im.hc = img;
    
    disp(rivmapsao_lourenco)
    
end

save('rivmapsao_lourenco')