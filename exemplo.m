RivMAP_path = 'C:\Users\renan\Unifesp\Iniciação Científica\ESTIMATIVA DA CONCENTRAÇÃO DE SEDIMENTOS EM SUSPENSÃO NA BACIA HIDROGRÁFICA DO RIO TAQUARIMS\RivMAP';
cd(RivMAP_path)
Ist = imread('exemplo.tif');
info = imfinfo('exemplo.tif');
img = logical(Ist);



riv_nosso.meta.year = 2020
riv_nosso.meta.exit_sides = 'NS';
riv_nosso.meta.Wn = 30;

riv_nosso.im = img;