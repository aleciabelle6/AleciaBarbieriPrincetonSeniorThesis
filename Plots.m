classdef Plots
    methods(Static)
        function subplot4againstAlpha(alpha, array1, array2, array3, array4,...
            YT1, YT2, YT3, YT4, mainTitle)
            figure;
            % Subplot 1: Panel (just wing)
            subplot(2,2,1);
            plot(alpha, array1, 'r-o', 'LineWidth', 1.5);
            title(YT1(2));
            xlabel('\alpha (deg)');
            ylabel(YT1(1));
            grid on;
            
            % Subplot 2: VLM (just wing)
            subplot(2,2,2);
            plot(alpha, array2, 'b-s', 'LineWidth', 1.5);
            title(YT2(2));
            xlabel('\alpha (deg)');
            ylabel(YT2(1));
            grid on;
            
            % Subplot 3: Panel (wing + human)
            subplot(2,2,3);
            plot(alpha, array3, 'm-^', 'LineWidth', 1.5);
            title(YT3(2));
            xlabel('\alpha (deg)');
            ylabel(YT3(1));
            grid on;
            
            % Subplot 4: Inferred VLM (wing + human)
            subplot(2,2,4);
            plot(alpha, array4, 'k-d', 'LineWidth', 1.5);
            title(YT4(2));
            xlabel('\alpha (deg)');
            ylabel(YT4(1));
            grid on;
        
            sgtitle(mainTitle)
        end
    end
end
