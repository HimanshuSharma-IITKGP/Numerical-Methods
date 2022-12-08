

function [] = timelinePlot(mat, x_arr)
    r = size(mat, 1);

    for i=1:r
        plot(x_arr, mat(r, :), linewidth=2);
        pause(0.2);
    end
end