clear
close all

% NPB7-1気球の気球重量 (皮膜重量+網重量)
W = 61.4178;
filename = 'NPB7';

% 網交点間隔
l0 = 0.102;

% 皮膜面密度と網線密度
rho_E = 9.20 * 10^-3;
rho_C = 3.64 * 10^-4;

% α=[0.1, sqrt(2/3)] の範囲で容積変化を調べる
a1 = linspace(0.1, 0.7, 200);
a2 = linspace(0.7, 0.8, 200);
a3 = linspace(0.8, sqrt(2/3), 3000);
% αの点列
% 重複を避けるため, a1(1:end-1), a2(end-1)とする.
% a3の最終成分を含めると発散するので a3(1:end-1)とする.
alphas = [a1(1:end-1), a2(1:end-1), a3(1:end-1)];

for i = 1:numel(alphas)
    alpha = alphas(i);
    Ws(i) = Wbar(alpha, l0);
    Rs(i) = sqrt(W/Ws(i));
    Ls(i) = Rs(i)   * integral(@(x)dLbar_dxbar(x,alpha),0,1);
    Ss(i) = Rs(i)^2 * integral(@(x)dSufbar_dxbar(x,alpha),0,1);
    Vs(i) = Rs(i)^3 * integral(@(x)dVolbar_dxbar(x,alpha),0,1);
    n = 2*pi*Rs(i) / (alpha*l0);
end

% % integral関数では a3(end-1)以降の積分が発散して解が求まらない.
% シンボリック計算を利用してがんばる.
a4 = linspace(a3(end-1), sqrt(2/3), 100);
a5 = linspace(a4(end-1), sqrt(2/3), 100);
% a4(91)に対しては解が求まらない.
% a5は解が得られる点から何点かをピックアップしている.
% 以上は自宅のHPで計算した場合であり, CPUによって結果が異なる可能性がある.
a45 = [a4([2:90,92:99]), a5([20,30,40,50,60,84,96])];

syms xbar

RelTol = 1e-3;

for i = 1:numel(a45)
    disp(i);
    a = a45(i);

    % 無次元化表面積の関数
    numer = 4*pi* sqrt(1-a^2) * xbar;
    delim = sqrt ( (1-a^2) - xbar.^4.*(1-a^2*xbar.^2) );
    funS  = numer ./ delim;

    % 無次元化全長の関数
    delim = sqrt ( (1-a^2) - xbar.^4.*(1-a^2*xbar.^2) ) .* sqrt (1-a^2.*xbar.^2);
    funL  = 2*sqrt(1-a^2) ./ delim;

    Sbar(i) = vpaintegral(funS,0,1,'RelTol',RelTol, 'MaxFunctionCalls', Inf);
    Lbar(i) = vpaintegral(funL,0,1,'RelTol',RelTol, 'MaxFunctionCalls', Inf);
    Wbar(i) = rho_E*Sbar(i) + (2*pi*rho_C/l0)*Lbar(i)/a;
    Rs2(i)  = sqrt(W/Wbar(i));
    Sbar(i) = Rs2(i)^2*Sbar(i);
    Lbar(i) = Rs2(i)*Lbar(i);


    % 複素数のゴミとりとクラスの実数化
    Sbar(i) = double( real( Sbar(i) ) );
    Lbar(i) = double( real( Lbar(i) ) );
    Wbar(i) = double( real( Wbar(i) ) );
    Rs2(i)  = double( real( Rs2(i) ) );

    % 無次元化容積の関数
    numer = 2*pi* xbar.^4 .* sqrt(1-a^2*xbar.^2);
    delim = sqrt ( (1-a^2) - xbar.^4.*(1-a^2*xbar.^2) );
    funV = numer ./ delim;
    Vbar = vpaintegral(funV,0,1,'RelTol',RelTol, 'MaxFunctionCalls', Inf);
    Vs2(i) = Rs2(i)^3*Vbar;
    % ゴミ取りと実数化
    Vs2(i) = double ( real ( Vs2(i) ) );
end

% 上記のループ終了後, なぜかSbar等が実数化されていない(sym型となっている).
% 代入式の左辺がsym型のため, sym型に型変換しているものと思われる.
% 以下で確実に実数化する
Sbar = double (Sbar);
Lbar = double (Lbar);
Vs2  = double (Vs2);
Rs2  = double (Rs2);

% 点列を作り直す
alphas = [alphas, a45];
Ss = [Ss, Sbar];
Ls = [Ls, Lbar];
Rs = [Rs, Rs2];
Vs = [Vs, Vs2];

% 高度と密度の対応表をつくっておく
% (高度60kmまで)
devs = 50000;
hkms = linspace(0,60,devs);
[Ts,as,ps,rhos] = atmoscoesa(1000*hkms);

% 網線軸力の計算
for i = 1:numel(alphas)
    % 気球重量の検算
    alpha   = alphas(i);
    n(i)    = (2*pi*Rs(i)) / (alpha*l0);
    Wes(i)  = rho_E * Ss(i);
    Wcs(i)  = rho_C * n(i) * Ls(i);
    Wall(i) = Wes(i) + Wcs(i);

    % 気球レベルフライト時の大気密度
    % ペイロード重量は気球構造重量と同等と仮定
    rho_air(i) = 2*Wall(i) / Vs(i);

    % レベルフライト時の高度/大気圧を求める
    ind   = min(find(rho_air(i) > rhos));
    hs(i) = hkms(ind);
    p     = ps(ind);

    % 必要差圧(大気圧の20%とする)
    dps = 0.2 * p;

    % 軸力
    c     = sqrt(1-alpha^2);  % c: cos(alpha)
    Ns(i) = pi*Rs(i)^2*dps / (n(i)*c);
end

figure
plot(alphas, Vs, '-o');
savefig([filename, '_V.fig']);

figure
plot(alphas, hs);
savefig([filename, '_H.fig']);

figure
plot(alphas, Ns);
savefig([filename, '_N.fig']);

figure
plot(alphas, Wall);

figure
plot(alphas, Rs);
