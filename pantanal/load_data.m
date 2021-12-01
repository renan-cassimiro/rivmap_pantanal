RivMAP_path = 'C:\Users\renan\Unifesp\Iniciação Científica\ESTIMATIVA DA CONCENTRAÇÃO DE SEDIMENTOS EM SUSPENSÃO NA BACIA HIDROGRÁFICA DO RIO TAQUARIMS\RivMAP\pantanal';
cd(RivMAP_path)

from = 1990;
to = 1996;


for year = from:to
    index = year-from + 1;

    img_name= strcat('channel_images/Active_channel_binary_mask_', string(year), '.tif');
    img_number = imread(img_name);
    img = logical(img_number);
    
    riv_pantanal(index).meta.year = year
    riv_pantanal(index).meta.exit_sides = 'NS';
    riv_pantanal(index).meta.Wn = 30;
    
    riv_pantanal(index).im = struct();
    riv_pantanal(index).im.st = img;
    riv_pantanal(index).im.hc = img;
    
    disp(riv_pantanal)
    
end

save('riv_pantanal')