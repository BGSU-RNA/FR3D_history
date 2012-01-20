% zIsoDiscrepancy(M1,M2,N1,N2) calculates the isodiscrepancy between the pairs represented by nucleotides M1,M2 and N1,N2

function [d] = zIsoDiscrepancy(M1,M2,N1,N2)

    R1 = M2.Rot' * M1.Rot;
    R2 = N1.Rot' * N2.Rot;

    ang = zAngleOfRotation(R1 * R2);

    t1 = (M2.Sugar(1,:) - M1.Sugar(1,:))*M1.Rot - (N2.Sugar(1,:) - N1.Sugar(1,:))*N1.Rot;
    t2 = (M1.Sugar(1,:) - M2.Sugar(1,:))*M2.Rot - (N1.Sugar(1,:) - N2.Sugar(1,:))*N2.Rot;

    cp = abs(norm(M2.Sugar(1,:) - M1.Sugar(1,:)) ...
           - norm(N2.Sugar(1,:) - N1.Sugar(1,:)));

    % calculate isodiscrepancy

    d = sqrt((4*ang)^2 + (t1*t1' + t2*t2')/2 + (3*cp)^2);
