clc
clearvars

load("var.mat")

% d ÄR EN DICTIONARY SOM INNEHÅLLER ALLA VÄRDEN PÅ KONSTANTER SOM ANVÄNDS FÖR ATT LÖSA LABB 1

function a = h(t, d)
% Genererar funktionen h'(t)

%Laddar konstanter
    H_1 = 0.22;
    H_2 = 0.4;
    v = d{"v"};
    L = d{"L"};
    if t <= L/v
        a = H_1/2 * (1 - cos(2*pi*v*t/L));
    elseif (7*L/v <= t) && (t <= 8*L/v)
        a = H_2/2 * (1 - cos(2*pi*v/L*(t - 7*L/v)));
    else
        a = 0;
    end
end

function a = h_prim(t, d)
    % Genererar funktionen h'(t)

%Laddar konstanter
    H_1 = 0.22;
    H_2 = 0.4;
    v = d{"v"};
    L = d{"L"};
    if t <= L/v
        a = H_1*pi*v/L * sin(2*pi*v*t/L);
    elseif (7*L/v <= t) && (t <= 8*L/v)
        a = H_2*pi*v/L * sin(2*pi*v/L*(t - 7*L/v));
    else
        a = 0;
    end
end

function [t_lista, eulerlista] = Euler(steg, steglangd, startvektor, d)
% Framåt Euler

%Laddar in konstanter
    k_2 = d{"k_2"};
    c_2 = d{"c_2"};
    M = d{"M"};
    A = d{"A"};
    t = 0;
    eulerlista = startvektor;
    t_lista = (0:steglangd:steglangd * steg)';
    for j = 1:steg
        F = [0; k_2*h(t, d) + c_2*h_prim(t, d)];
        g = [0; 0; M\F];
        w_innan = eulerlista(:,j);
        w_efter = w_innan + steglangd*(A*w_innan + g);
        eulerlista(:, j+1) = w_efter;
        t = t + steglangd;
    end
end

function dydt = funktion1(t, Y, d)
% System av differentialekvationer till ODE45
    A = d{"A"};
    M = d{"M"};
    k_2 = d{"k_2"};
    c_2 = d{"c_2"};

    %Detta matas in i ODE45-funktionen
    dydt = A * [Y(1), Y(2), Y(3), Y(4)]' + [0; 0; M\[0; k_2*h(t, d) + c_2*h_prim(t, d)]];
end

function delta_t = var_tidssteg(t)
% Hjälpfunktion för att finna steglängd från tidlista. Tidslistan genereras
% av ODE-45-funktionen
   delta_t = zeros(length(t)-1, 1);
   for i = 1:(length(t)-1)
        delta = t(i+1) - t(i);
        delta_t(i) = delta;
   end
end

function a = stab_villkor(d)
% Hittar minsta tidssteg för en diagonaliserbar matris A för vilket
% systemet i frågan är stabilt
    A = d{"A"};
    egen = eig(A);
    size = length(egen);
    a = Inf;
    for i = 1:size
        ev = egen(i, 1);
        temp = (-2)*real(ev)/(real(ev)^2 + imag(ev)^2);
        if temp < a
            a = temp;
        end
    end
end

function d = nya_d(varden)
% Ändrar på värden i d
% varden är en cell-array som innehåller värdena som motsvarar konstanterna
% som namnges i namn
    namn = ["m_1", "m_2", "k_1", "k_2",  "c_1", "c_2", "v", "H", "L"];
    d = dictionary(namn, varden);

    M = [d{"m_1"}, 0; 0, d{"m_2"}]; 
    K = [d{"k_1"}, (-1)*d{"k_1"}; (-1)*d{"k_1"}, d{"k_1"} + d{"k_2"}];
    C = [d{"c_1"}, (-1)*d{"c_1"}; (-1)*d{"c_1"}, d{"c_1"} + d{"c_2"}];

    M_K = (-1)*M\K;
    M_C = (-1)*M\C;

    A = [0, 0, 1, 0; 
         0, 0, 0, 1; 
         M_K, M_C];

    d = insert(d, "M", {M});
    d = insert(d, "A", {A});
end

function a = F(t, d)
% F-funktion i högerledet från uppgiftsbeskrivningen. a är funktionsvärdet
% vid en given tidpunkt. 
% Används endast i implicita Heunmetoden
    k_2 = d{"k_2"};
    c_2 = d{"c_2"};
    a = [0; k_2*h(t, d) + c_2*h_prim(t , d)];
end

function [t_lista, heunlista] = imp_heun(steg, steglangd, startvektor, d)
% Implicita Heuns metod / implicita trapetsmetoden
    A = d{"A"};
    M = d{"M"};
    heunlista = startvektor;
    t_lista = (0:steglangd:steglangd * steg)';
    for j = 1:steg
        % B_1 och B_2 är konstanta matriser från härledningen av implicita
        % trapetsmetoden.
        t = t_lista(j);
        B_1 = (2/steglangd)*eye(4) - A; % B_1 = (2/h*I - A)
        B_2 = (2/steglangd)*eye(4) + A; % B_2 = (2/h*I - A) 
        g_1 = [0; 0; M\F(t, d)]; % g_n från uppgiften
        g_2 = [0; 0; M\F(t + steglangd, d)]; % g_n+1 från uppgiften
        w_innan = heunlista(:,j);
        w_j = B_1\(B_2*w_innan + (g_1 + g_2)) ;
        heunlista(:,j+1) = w_j; 
    end
end

function maxvarden = flera_numeriska_losningar(d, skalfaktorer, metod, h_0, langdskala)
% Kör alla metoder med olika värden på tidssteget som motsvarar
% alfa*delta_tmax, där delta_tmax är stabilitetsvillkoret för framåt Euler
% h_0 är stabilitetsvillkoret
% langdskala är intervallet vi plottar i
    size = length(skalfaktorer);
    maxvarden = zeros(size, 1);

    for i = 1:size
        a = skalfaktorer(1, i);
        ftemp = figure(name="alfa="+a);
        
        if metod == "E"
            [t_temp, losningar_temp] = Euler(round(langdskala/(a*h_0)), a*h_0, [0, 0, 0, 0]', d);
        elseif metod == "T"
            [t_temp, losningar_temp] = imp_heun(round(langdskala/(a*h_0)), a*h_0, [0, 0, 0, 0]', d);
        end
        figure(ftemp);
        plot(t_temp, losningar_temp(2, :))
    end
end

function N = Newtonsmethod(k_1_ref, k_2_ref, n)
% Newtons metod i flera variabler
% n är ett "back-up" - stoppvillkor: Om inte toleransen nås inom rimlig tid
% så stoppas programmet istället efter n iterationer
    gissning = [k_1_ref;k_2_ref];
    tol = 10^-6;

    for i = 1:n
        F = transfer_functions(gissning(1), gissning(2), 0.81, 1);
        J = Jacobian_transfer_functions(gissning(1), gissning(2));
        temp = gissning;
        gissning = gissning - J\F; 
        if abs(gissning(1)-temp(1)) < tol && abs(gissning(2)-temp(2)) < tol
            break
        end
    end
    N = gissning; % N är den slutgiltiga gissningen på k1 och k2
end

function list = konv(error_list, m) % m = konvergensordning
% Funktion från labb 2 för att hitta konvergensordningen givet en lista med
% fel.
    n = length(error_list);
    list = zeros(1, n-1);
    for i = 1:(n-1)
        list(i) = error_list(i+1)/(error_list(i)^m);
    end
end

%% Följande funktioner löser uppgifterna:

function dummy = uppgift2(d)
% Uppgift 2
    f_1 = figure(name="ODE45");
    f_2 = figure(name="Tidstegsplot");
    f_3 = figure(name="Euler vs. ODE45");

    % a)
    disp("2. a)")
    
    % Parametrar för ODE45:s inställningar
    options = odeset("RelTol", 10^-6);
    [t, y] = ode45(@(t, Y) funktion1(t, Y, d), [0, 1], [0,0,0,0]', options);
    
    % Maxvärden för utslagen i z1 och z2
    max_z1 = max(abs(y(:, 1)));
    max_z2 = max(abs(y(:, 2)));

    disp(['maxvärdet för z_1: ', num2str(max_z1), ' m'])
    disp(['maxvärdet för z_2: ', num2str(max_z2), ' m'])

    figure(f_1)
    plot(t, y(:,1), "b-")
    hold on
    plot(t, y(:,2), "r-")
    legend("z_1", "z_2")

    % b)

    % Parametrar för ODE45:s inställningar
    steg_options = odeset("RelTol", 10^-6, "Refine", 1);
    % t_steg i output =/= t_steg i funktionen i ODE45:s input
    [t_steg, y_dummy] = ode45(@(t_steg, Y) funktion1(t_steg, Y, d), [0, 1], [0,0,0,0]', steg_options);

    delta_t = var_tidssteg(t_steg);

    figure(f_2);
    plot(1:length(delta_t), delta_t, "ro")
    legend("Längd på tidssteg (m)")

    % c)

    % Löser med Euler för olika h
    [t_lista_1, eulerlista_1] = Euler(200, 5*10^-3, [0, 0, 0, 0]', d);
    [t_lista_2, eulerlista_2] = Euler(10000, 10^-4, [0, 0, 0, 0]', d);

    figure(f_3);
    plot(t_lista_1, eulerlista_1(2, :), "r-");
    hold on 
    plot(t_lista_2, eulerlista_2(1, :), "b-");
    plot(t, y(:,2), "g-")

    legend("Euler 5*10^-3", "Euler 10^-4", "ODE45")
    max(abs(eulerlista_2(1, :)))
end

%%

function dummy = uppgift3(d)

    % b)
    disp("3. b)")

    delta_tmax = stab_villkor(d);
    disp(['Maximala tidsteget för stabilitet givet tabellvärdena är: ', num2str(delta_tmax)])

    uppfyller_stabilitetsvillkor = delta_tmax > 5*10^-3

    % c)
    skalfaktorer = [0.9, 1, 1.1, 1.5]; % Lista med "alfan"
    
    flera_numeriska_losningar(d, skalfaktorer, "E", delta_tmax, 2);

    % d)
    disp("3. d)")

    % k2 görs 100 gånger styvare
    d2 = nya_d({475,    53,    5400,  135000*100, 310,   1200,  65/3.6, 0.24, 1});
    f1 = figure(name="Hårdare fjädring");
    
    % Nytt stabilitetsvillkor för system med styvare fjäder
    delta_tmax_2 = stab_villkor(d2)

    [t, eulerlista] = Euler(10000, 0.1*delta_tmax_2, [0, 0, 0, 0]', d2);

    figure(f1);
    plot(t, eulerlista(1, :), "b-")
    hold on
    plot(t, eulerlista(2, :), "r-")
    legend("z_1", "z_2")
end

%%

function dummy = uppgift4(d)
    %a)
    d = nya_d({475,    53,    5400,  135000*100, 310,   1200,  65/3.6, 0.24, 1});
    skalfaktorer1 = [1, 10, 100];
    delta_tmax = stab_villkor(d);

    flera_numeriska_losningar(d, skalfaktorer1, "T", 10^-3, 0.065);

    % b)
    disp("4. b)")
    delta_t_noll = 10^-3;
    
    skalfaktorer2 = [1, 1/2, 1/4, 1/8];
    % Lista av maximala felen för z2 för respektive alfa i skalfaktorer2
    maxvarden = zeros(1, length(skalfaktorer2));
    langdskala = 0.05;

    for i = 1:length(skalfaktorer2)
        a = skalfaktorer2(i);
        options = odeset("RelTol", 10^-9, "AbsTol",10^-9);
        [t_temp, losningar_temp] = imp_heun(round(langdskala/(a*delta_t_noll)), a*delta_t_noll, [0, 0, 0, 0]', d);
        [t, y] = ode45(@(t, Y) funktion1(t, Y, d), t_temp, [0,0,0,0]', options);

        absvarde_temp = abs(y(:, 2) - losningar_temp(2, :)');
        maxvarde_temp = max(absvarde_temp);
        maxvarden(i) = maxvarde_temp;
    end

    for i = 1:(length(skalfaktorer2)-1)
        % Beräknar approximativ noggrannhetsordning för Heuns metod med
        % respektive skalfaktor alfa
        p_temp = log(maxvarden(i)/maxvarden(i+1))/log(2);
        disp("Beräknad noggrannhetsordning utifrån alfa1 = " + skalfaktorer2(i) + ...
            " och alfa2 = " + skalfaktorer2(i+1) + ": " + p_temp)
    end
    
    format longG
    disp("Maximala absbeloppen för skalfaktorer:");
    disp(maxvarden)
    % Eftersom globala felen blir begränsade när vi delar med h^2 vet vi
    % att noggrannhetsordningen är 2, vilket stämmer överens med teorin för
    % trapetsmetodens noggrannhetsordning
end

%%

function dummy = uppgift5(d)
    % a)
    disp("5. a)")
    k_1_ref = d{"k_1"};
    k_2_ref = d{"k_2"};
    N = Newtonsmethod(k_1_ref, k_2_ref, 100);
    
    k_1 = N(1);
    k_2 = N(2);

    disp("De nya värdena på k_1 och k_2 blir:")
    disp("k_1 = " + k_1)
    disp("k_2 = " + k_2)

    % b)
    N = Newtonsmethod(d{"k_1"}, d{"k_2"}, 1000);
    a = zeros(100, 2);
    f1 = figure(name="felverifiering");
    
    % OBS!!! Dumheter börjar här
    for i = 1:100
        e_1 = Newtonsmethod(d{"k_1"}, d{"k_2"}, i);
        e_2 = Newtonsmethod(d{"k_1"}, d{"k_2"}, i+1);
        fel_1 = norm(e_1 - N);
        fel_2 = norm(e_2 - N);
        a(i, 1) = fel_1;
        a(i, 2) = fel_2;
    end

    %Kontrollerar att felet är kvadratiskt
    figure(f1)
   
    x = -10:10;
    %Approoximerar log(M) med 9
    y = 2*x - 9;
    hold on
    plot(log(a(:,1)), log(a(:,2)), "g*")
    plot(x, y, "r")

    legend("logaritmerade kvadratfelsdivisionen", "y = 2x - ~log(M)")
    % OBS!!! Dumheterna upphör

    % c)

    % Nya värden på k1 och k2 sätts in i d2
    d2 = nya_d({475, 53, k_1, k_2, 310, 1200, 65/3.6, 0.24, 1});
    f2 = figure(name="z1 ref vs. z1 optimal");

    % Ser till att tidssteget är inom konvergensradien för A i d och d2
    delta_tmax1 = stab_villkor(d);
    delta_tmax2 = stab_villkor(d2);
    delta_tmax = min([delta_tmax1, delta_tmax2]);
    
    [t_lista1, eulerlista1] = Euler(round(2/(delta_tmax*0.1)), delta_tmax*0.1, [0, 0, 0, 0]', d);
    [t_lista2, eulerlista2] = Euler(round(2/(delta_tmax*0.1)), delta_tmax*0.1, [0, 0, 0, 0]', d2);
    
    % Plottar z1 med ursprungliga k1 och k2 samt med optimerade k1 och k2
    figure(f2)
    plot(t_lista1, eulerlista1(1, :), "r")
    hold on
    plot(t_lista2, eulerlista2(1, :), "b")

    legend("k1 ref", "k1 optimal")
end

uppgiftsval = input("Välj uppgift. 0 för att köra alla");

if uppgiftsval == 2
    uppgift2(d)
elseif uppgiftsval == 3
    uppgift3(d)
elseif uppgiftsval == 4
    uppgift4(d)
elseif uppgiftsval == 5
    uppgift5(d)
elseif uppgiftsval == 0
    uppgift2(d)
    uppgift3(d)
    uppgift4(d)
    uppgift5(d)
end