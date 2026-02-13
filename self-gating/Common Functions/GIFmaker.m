function GIFmaker(figHandle, gifFileName, delayTime)
    % Ensure the figure is rendered on-screen
    drawnow;

    % Capture the current figure as an image
    frame   = getframe(figHandle);
    img     = frame2im(frame);

    % Convert the image to indexed format
    [imind, cm] = rgb2ind(img, 256);

    % Check if the GIF file already exists
    if exist(gifFileName, 'file')
        % Append the image to the existing GIF file
        imwrite(imind, cm, gifFileName, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    else
        % Create a new GIF file with the current image as the first frame
        imwrite(imind, cm, gifFileName, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    end
end