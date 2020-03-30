%{
    Ricardo Lopez R. A01066515                                                  27/03/2020    
    Héctor Irwin M. A01746543
    Donaldo Garrido A01275416
    
    Espectro de Energía de la Molécula Mg2 - (Tiempo aproximado de ejecución: 287.052427 s)

       - El siguiete programa utiliza distintos métodos númericos para el
        cálculo de espectro de energía de la molecula Mg2.
        
        - Lo resultablos de de las distintas energpias se imprimen en la
        consola al terminar el programa para l = {0,1,2}. Además se
        despliega la gráfica del potencial efectivo con sus respectivas
        energías donde l = 1 y l =2
        
%}

%Funcion principal
function energySpectrumMg2()

    format long
    close all
    clear all
    
    %Calculo de energias cuando l = {0,1,2}
    tic
    for i = 1 : 3
      eigenE(:,i) = mainenergySpectrumMg2(i) ;    
    end
    
    %Funcion para graficar 
    ploting(eigenE);
    toc
end 

function eigenE =  mainenergySpectrumMg2(num)

  %Declaracion de constantes y potencial
    a = 4.96 * 10^7;
    b = 624;
    m1 = 24*1836 + 12 ;
    m2 =24*1836 + 12 ;
    r = linspace(6,19,100000); 
    M = (m1*m2)/(m1 + m2);  
    l = num - 1;
    Vr = (a./r.^12) - (b./r.^6) + (l^2 + l + 1/4)*(1./(2*M* r.^2));    
    
    intervalos = intervals(Vr,length(Vr),r);
    
    disp("CALCULANDO ENERGIAS para l = "+l+"........")
    index1 = find(r == intervalos(2));
    index2 = find( r > 18.9 & r < 19,1);
    
    count = 1;
    %Creacion de un vector E de energias en el rango de interes
    for i = index1 : index2
        E(count) = Vr(i);
        count = count + 1;
    end        
   
    roots = zeros(length(E)-1,2);   
    
    for j = 1 : length(E) - 1
        
        f = E(j) - Vr;

    %--------Busqueda de intervalos f = E - V(r)-----------%
        inter = intervals(f,length(f),r);
        dimensions = size(inter);      

        if size(inter,1) < 2        
            root1 = r(f == 0);
            roots(j,1) = root1;  
        end
        
    %------------Busqueda de todas las raices por  Bisecccion-----%
        for i = 1: dimensions(1)
            roots(j,2) = bisection(inter(i,1),inter(i,2),E(j),l);
        end
    end   
    
    I = zeros(1,size(roots,1));
    
    for i = 1 : size(roots,1)
        
        points = 1000;
        r1 = roots(i,1);
        r2 = roots(i,2);

        R = linspace(r1,r2,points);
        Vr2 = (a./R.^12) - (b./R.^6) + (l^2 + l + 1/4)*(1./(2*M* R.^2)) ;
        F = sqrt( E(i) - Vr2);    

        R2 = linspace(r1,r2,points*2);
        Vr3 = (a./R2.^12) - (b./R2.^6) + (l^2 + l + 1/4)*(1./(2*M* R2.^2)) ;
        F2 = sqrt( E(i) - Vr3);   

        I1 = trapz(R,F);
        I2 = trapz(R2,F2);
        I(1,i) = (sqrt(2*M)/pi) * rombergIntegration(I1,I2);        %Calculo de las integrales para las distintas energias
    end
    
    EBK = zeros(1,size(I,2));
    
    for j = 1 : 11
        
        n = j -1;
        contador = 1;
        disp("----------------- N = "+n+" ||  l = "+l+" ----------------------");      
        
        for i = 1 : size(I,2)   
            
            EBK(i) = I(i) - (n + 0.5);    % f(E) = EBK = 0
            %Buscando cual energia aproxima mas la funcion a 0 (Primer filtro)
            if EBK(i) >= -0.009 && EBK(i) <= 0.009
                
                eigenEnergy(contador) = E(i);
                ebkList(contador)  = abs(EBK(i));
                contador = contador + 1;      
                
            end      
            
        end
        
        minIndex = ebkList == min(ebkList);
        energy(j) = eigenEnergy(minIndex);        
        disp("E_n = "+eigenEnergy(minIndex));    
        eigenEnergy = [];
        ebkList = [];
        
    end    
    eigenE = energy;
    energy = [];
end

%-----------Busqueda de intervalos----------------%
function [intervalos] =  intervals(f,length,x)

    count = 1;
    bool = false;
    for i = 1:length-1

        if f(i)*f(i+1) < 0 
            bool = true;
            xUpper = x(i+1);
            xLower = x(i);
            interval(count,1) = xLower;
            interval(count,2) = xUpper;
            count = count + 1;
        end
    end 
    
    if bool == true
        intervalos = interval;
        
    else 
        intervalos = [];
    end
end


%-----------Metodo de Biseccion--------------%
function [root] = bisection(xUpper,xLower,E,L)  

    a = 4.96 * 10^7;
    b = 624;
    m1 = 24*1836 + 12 ;
    m2 =24*1836 + 12 ;
    M = (m1*m2)/(m1 + m2);        
    l = L;
    
    error = 100;
    xM =  (xLower + xUpper)/2;
    
    while error > 0.000001 
        
        newxUpper = xM;
        r = [xLower,newxUpper];
        
        %F(x) declaration
        f = E - ( (a./r.^12) - (b./r.^6) + (l^2 + l + 1/4)*(1./(2*M* r.^2)) );
        
        if f(1)*f(2) < 0             
            xUpper = xM;
            newxM = (xLower+xUpper)/2; 
        end
        
        if f(1)*f(2) > 0          
            xLower = xM;
            newxM = (xLower+xUpper)/2 ;
        end
        
        error = ((abs(newxM - xM))*100)/(newxM);
        xM = newxM;
        
        if f(1)*f(2) == 0
            disp("La raiz es xM");
            error = .00001;
        end 
        
    end    
    root = xM;
end

%--------Integracion De Romberg (Trapz con Extrapolacion de Richardson)------%
function [integral] = rombergIntegration(I1,I2)
    I = (4/3)*(I2) - (1/3)*(I1);    
    integral = I;
end

%Funcion para  graficar 
function ploting(energies)

    energy = energies;     
    a = 4.96 * 10^7;
    b = 624;
    m1 = 24*1836 + 12 ;
    m2 =24*1836 + 12 ;
    r = linspace(0,35,1000); 
    M = (m1*m2)/(m1 + m2);  
    %graficado para l = 1 y l = 2
    for i = 1 : 2
        
        l = i;        
        Vr = (a./r.^12) - (b./r.^6) + (l^2 + l + 1/4)*(1./(2*M* r.^2));  
        
        figure(i)
        plot(r,Vr,'b');
        title("Potencial efectivo V(r) , donde l = "+i)
        xlabel('$\textbf{r}$','interpreter','latex')
        ylabel('$\textbf{V(r)}$','interpreter','latex')    
        xlim([0 15])
        ylim([-2e-3,2e-3]);
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin'; 
        hold on    
        
        for j = 1 : size(energy,1)    
           hi = plot(r,ones(size(r))*energy(j,i+1)) ;
           h(j) = hi(1);
           if j == 1
            hold on 
           end
           LegendsStrings{j} = ['E_ ',int2str(j-1),' = ',num2str(energy(j,i+1))];
        end
        legend(h,LegendsStrings, 'Interpreter', 'none');
        LegendsStrings = [];
    end
end