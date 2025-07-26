classdef AeroCoeffs
    methods(Static)
        function [C_L, C_D, C_di, C_d0, L_D, eff] = AeroCoeffsFromAlpha(alphaTarget, alphaArr, C_LArr, C_DArr, C_diArr, C_d0Arr, L_DArr, effArr)
            [C_L] = AeroCoeffs.YfromAlpha(alphaArr, C_LArr, alphaTarget);
            [C_D] = AeroCoeffs.YfromAlpha(alphaArr, C_DArr, alphaTarget);
            [C_di] = AeroCoeffs.YfromAlpha(alphaArr, C_diArr, alphaTarget);
            [C_d0] = AeroCoeffs.YfromAlpha(alphaArr, C_d0Arr, alphaTarget);
            [L_D] = AeroCoeffs.YfromAlpha(alphaArr, L_DArr, alphaTarget);
            [eff] = AeroCoeffs.YfromAlpha(alphaArr, effArr, alphaTarget);
        end
        
        function [targetY] = YfromAlpha(alphaArr, yArr, alphaTarget)
            pY = polyfit(alphaArr, yArr, 3);
            targetY = polyval(pY, alphaTarget);
        end
    end
end