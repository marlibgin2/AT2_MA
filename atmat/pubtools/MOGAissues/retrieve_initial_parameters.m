function X0 = retrieve_initial_parameters(rrr, tipo)

if nargin < 2
    tipo = 'B';
end

switch tipo
    case 'B'
        q1i = findcells(rrr,'FamName','Q1');
        q2i = findcells(rrr,'FamName','Q2');
        q3i = findcells(rrr,'FamName','Q3');
        q4i = findcells(rrr,'FamName','Q4');
        q5i = findcells(rrr,'FamName','Q5');
        q6i = findcells(rrr,'FamName','Q6');

        r1i = findcells(rrr,'FamName','R1');
        r2i = findcells(rrr,'FamName','R2');

        s1i = findcells(rrr,'FamName','S1');
        s2i = findcells(rrr,'FamName','S2');
        s3i = findcells(rrr,'FamName','S3');
        s4i = findcells(rrr,'FamName','S4');
        s5i = findcells(rrr,'FamName','S5');

        o1i = findcells(rrr,'FamName','O1');
        o2i = findcells(rrr,'FamName','O2');
        o3i = findcells(rrr,'FamName','O3');

        dipi   = [findcells(rrr,'FamName','D2') findcells(rrr,'FamName','D3')];
        dipmi  = findcells(rrr,'FamName','D1');

        %     for i=1:length(dipi)
        %         wdip(i) = rrr{dipi(i)}.PolynomB(2)/rrr{dipi(7)}.PolynomB(2);
        %     end
        %     for i=1:length(dipmi)
        %         wdipm(i) = rrr{dipmi(i)}.PolynomB(2)/rrr{dipmi(7)}.PolynomB(2);
        %     end

        VARi = {q1i; q2i; q3i; q4i; q5i; q6i}; nVAR(1) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j) = rrr{VARi{j}(1)}.PolynomB(2);
        end

        VARi = {r1i; r2i}; nVAR(2) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(2);
        end

        clear VARi;
        VARi = {dipi; dipmi}; nVAR(3) = numel(VARi);
        %     W    = {wdip; wdipm};
        for j = 1:2
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(7)}.PolynomB(2);
        end

        clear VARi;
        VARi = {s1i; s2i; s3i; s4i; s5i}; nVAR(4) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(3);
        end

        clear VARi;
        VARi  = {o1i; o2i; o3i}; nVAR(5) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(4);
        end

    case 'B+'
        q1i = findcells(rrr,'FamName','Q1');
        q2i = findcells(rrr,'FamName','Q2');
        q3i = findcells(rrr,'FamName','Q3');
        q4i = findcells(rrr,'FamName','Q4');
        q5i = findcells(rrr,'FamName','Q5');
        q6i = findcells(rrr,'FamName','Q6');

        r1i = findcells(rrr,'FamName','R1');
        r2i = findcells(rrr,'FamName','R2');

        s1i = findcells(rrr,'FamName','S1');
        s2i = findcells(rrr,'FamName','S2');
        s3i = findcells(rrr,'FamName','S3');
        s4i = findcells(rrr,'FamName','S4');
        s5i = findcells(rrr,'FamName','S5');

        o1i = findcells(rrr,'FamName','O1');
        o2i = findcells(rrr,'FamName','O2');
        o3i = findcells(rrr,'FamName','O3');
               
        dipi   = findcells(rrr,'FamName','D2');
        dip3i  = findcells(rrr,'FamName','D3');
        dipmi  = findcells(rrr,'FamName','D1');

          
        VARi = {q1i; q2i; q3i; q4i; q5i; q6i}; nVAR(1) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j) = rrr{VARi{j}(1)}.PolynomB(2);
        end

        VARi = {r1i; r2i}; nVAR(2) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(2);
        end
      
        clear VARi;
        VARi = {dipi; dip3i; dipmi}; nVAR(3) = numel(VARi);
        %     W    = {wdip; wdipm};
        for j = 1:3
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(7)}.PolynomB(2);
        end

                clear VARi;
        VARi = {s1i; s2i; s3i; s4i; s5i}; nVAR(4) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(3);
        end

        clear VARi;
        VARi  = {o1i; o2i; o3i}; nVAR(5) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(4);
        end

       
    case 'B03'
        q1i = findcells(rrr,'FamName','Q1');
        q2i = findcells(rrr,'FamName','Q2');
        q3i = findcells(rrr,'FamName','Q3');
        q4i = findcells(rrr,'FamName','Q4');

        r1i = findcells(rrr,'FamName','R1');

        s1i = findcells(rrr,'FamName','S1');
        s2i = findcells(rrr,'FamName','S2');
        s3i = findcells(rrr,'FamName','S3');
        s4i = findcells(rrr,'FamName','S4');
        s5i = findcells(rrr,'FamName','S5');

        o1i = findcells(rrr,'FamName','O1');
        o2i = findcells(rrr,'FamName','O2');
        o3i = findcells(rrr,'FamName','O3');
               
        dipi   = findcells(rrr,'FamName','D2');
        dip3i  = findcells(rrr,'FamName','D3');
        dipmi  = findcells(rrr,'FamName','D1');

          
        VARi = {q1i; q2i; q3i; q4i}; nVAR(1) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j) = rrr{VARi{j}(1)}.PolynomB(2);
        end

        VARi = {r1i}; nVAR(2) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(2);
        end
      
        clear VARi;
        VARi = {dipi; dip3i; dipmi}; nVAR(3) = numel(VARi);
        %     W    = {wdip; wdipm};
        for j = 1:3
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(7)}.PolynomB(2);
        end

                
        clear VARi;
        VARi = {s1i; s2i; s3i; s4i; s5i}; nVAR(4) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(3);
        end

        clear VARi;
        VARi  = {o1i; o2i; o3i}; nVAR(5) = numel(VARi);
        for j = 1:numel(VARi)
            X0(j+sum(nVAR(1:end-1))) = rrr{VARi{j}(1)}.PolynomB(4);
        end

        
end



end