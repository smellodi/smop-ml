% Creates a graph to replace Bubble-plot from MOD_5-search

function utils()
    solution = [
        % copied from MOD_5-search
        10	42	3.793	9	43	1.41	-1	1	1.41	9	43	8.21;
        18	29	3.231	33	26	15.30	15	-3	2.24	16	28	12.48;
        30	20	2.9679	26	21	4.12	-4	1	1.00	30	21	8.50;
        39	8	2.7194	42	6	3.61	3	-2	2.00	39	6	21.70;
        11	12	5.2532	27	9	16.28	16	-3	5.00	6	12	12.18;
        37	38	4.2999	44	37	7.07	7	-1	3.61	39	41	13.87;
        3	26	2.4775	6	31	5.83	3	5	3.16	6	25	15.52;
        22	3	3.4022	18	4	4.12	-4	1	1.00	23	3	13.47;
        47	22	5.4886	40	24	7.28	-7	2	3.61	45	25	15.35;
        26	47	5.7019	26	50	3.00	0	3	1.00	26	46	10.11
    ];

    c1 = [1 0.6 0.3];
    c2 = [0.1 0.2 0.5];

    figure('Position',[100,100,400,400]);
    scatter(-1,[-1; -1],[],[c1; c2],"filled");
    plotEllipses(solution,4,5,c1);
    plotEllipses(solution,10,11,c2);

    xlabel('Isopropanol, sccm');
    ylabel('Ethanol, sccm');
    xlim([0,50]);
    ylim([0,50]);
    xticks(0:5:50);
    yticks(0:5:50);

    grid on;
    lgd = legend('Solution', 'Closest tested');
    lgd.Location = 'northoutside';
    lgd.NumColumns = 2;
    lgd.Position = [0.25 0.94 0.5 0.05];
    lgd.EdgeColor = [1 1 1];
    lgd.FontSize = 14;

    ax = gca;
    ax.FontSize = 12;
end

function plotEllipses(data,wi,hi,color)
    for jj = 1:size(data,1)
        line = data(jj,:);
        dx = abs(line(1) - line(wi));
        dy = abs(line(2) - line(hi));
        line(1) = line(1) - dx/2;
        line(2) = line(2) - dy/2;
        h = rectangle('Position',[line(1),line(2),dx,dy],'Curvature',[1,1]); 
        h.FaceColor = color;
        h.EdgeColor = color; 
        h.LineWidth = 1;
    end
end

