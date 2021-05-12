function [wynik] = proguj (obraz,prog)

[rx, ry] = size (obraz);

wynik = zeros (rx, ry);

for (x=1:rx)
    for (y=1:ry)
        if (obraz(x,y)>prog)
            wynik(x,y)=255;
        else
        wynik(x,y)=0;
        end
    end
end