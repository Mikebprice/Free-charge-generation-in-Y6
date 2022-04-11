function h = fig
ImageDPI=600; ImageSizeX=10.2; ImageSizeY=7;
ImageFontSize=10;
 h=figure('Units','centimeters','InnerPosition',[10 5 ImageSizeX ImageSizeY],...
          'PaperPosition',[10 5 ImageSizeX ImageSizeY])
map=0.75*colormap(cool(4));
end

