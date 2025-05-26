function txtHandle = TextZoomable(x,y,varargin)
% txtHandle = FixedSizeText(x,y,varargin)
%
% Version 2.0, 14 July 2015
% Adds text to a figure in the normal manner, except that this text
% grows/shrinks with figure scaling and zooming, unlike normal text that
% stays at a fixed font size during figure operations. Note it scales with
% figure height - for best scaling use 'axis equal' before setting up the
% text.
%
% All varargin{:} arguments will be passed directly on to the text
% function (text properties, etc.)
%
% (doesn't behave well with FontUnits = 'normalized')
%
% example:
%
% figure(1); clf;
% rectangle('Position', [0 0 1 1]);
% rectangle('Position', [.25 .25 .5 .5]);
% 
% th = TextZoomable(.5, .5, 'red', 'color', [1 0 0], 'Clipping', 'on');
% th2 = TextZoomable(.5, .1, 'blue', 'color', [0 0 1]);
%
%
% Ken Purchase, 4-25-2013, with many thanks to the Matlab User
% Community, including Matt J, Hoi Wong, Philip Caplan.
%
% Modified by Philip Caplan on 12/10/2014 to support R2014b
%
% Modified by Hoi Wong on 1/23/2014 so that the 'UserData' will not conflict
% with what is already set by gscatter() if previously called.
% create the text
%
% Modified by Ken, incorporating above changes and version compatibility
% for previous versions
%
% Note there is an issue on my computer that I haven't yet fixed -
% initially the font displays in an absurdly wrong size, but as soon as I
% zoom in or out, all is correct. I don't have time to dig into it now, but
% if you fix it, give me a comment on the TextZoomable page at the file
% exchange, and I can publish your fix. Many thanks.
%
%    
    
    txtHandle = text(x,y,varargin{:});
    % NOTE: The function can be disabled by returning here.
    
    % detect its size relative to the figure, and set up listeners to resize
    % it as the figure resizes, or axis limits are changed.
    hAx = gca;
    hFig = get(hAx,'Parent');
    
    fs = get(txtHandle, 'FontSize');
    ratios = fs * diff(get(hAx,'YLim')) / max(get(hFig,'Position') .* [0 0 0 1]);
    
    % append the handles and ratios to the user data - repeated calls will
    % add each block of text to the list
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      ud = getappdata(hAx, 'TextZoomable_UserData');
    else
      ud = get(hAx, 'UserData');
    end
    
    if isfield(ud, 'ratios')
      ud.ratios = [ud.ratios(:); ratios];
    else
      ud.ratios = ratios;
    end
    if isfield(ud, 'handles')
      ud.handles = [ud.handles(:); txtHandle];
    else
      ud.handles = txtHandle;
    end
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      setappdata(hAx, 'TextZoomable_UserData', ud);
    else
      set(hAx,'UserData', ud);
    end
    
    localSetupPositionListener(hFig,hAx);
    localSetupLimitListener(hAx);
end
  

%% Helper Functions
  function fs = getBestFontSize(imAxes)
    % Try to keep font size reasonable for text
    hFig = get(imAxes,'Parent');
    hFigFactor = max(get(hFig,'Position') .* [0 0 0 1]);
    axHeight = diff(get(imAxes,'YLim'));
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      ud = getappdata(imAxes, 'TextZoomable_UserData');
    else
      ud = get(imAxes,'UserData');  % stored in teh first user data.
    end
    
    fs = round(ud.ratios * hFigFactor / axHeight);
    fs = max(fs, 3);
  end
  
  
  function localSetupPositionListener(hFig,imAxes)
    % helper function to sets up listeners for resizing, so we can detect if
    % we would need to change the fontsize
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      PostPositionListener = addlistener(hFig,'SizeChanged',...
        @(o,e) localPostPositionListener(o,e,imAxes) );
    else
      PostPositionListener = handle.listener(hFig,'ResizeEvent',...
        {@localPostPositionListener,imAxes});
    end
    
    setappdata(hFig,'KenFigResizeListeners',PostPositionListener);
  end
  
  
  function localPostPositionListener(~,~,imAxes)
    % when called, rescale all fonts in image
    if( ~isvalid(imAxes) ) % The imAxes might be deleted, which getappdata() doesn't have an output and will throw an exception
      return;
    end
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      ud = getappdata(imAxes, 'TextZoomable_UserData');
    else
      ud = get(imAxes,'UserData');
    end
    
    fs = getBestFontSize(imAxes);
    for ii = 1:length(ud.handles)
      set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
    end
  end
  
  
  function localSetupLimitListener(imAxes)
    % helper function to sets up listeners for zooming, so we can detect if
    % we would need to change the fontsize
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      LimListener = addlistener(imAxes,{'XLim','YLim'},'PostSet',@localLimitListener);
    else
      hgp     = findpackage('hg');
      axesC   = findclass(hgp,'axes');
      LimListener = handle.listener(imAxes,[axesC.findprop('XLim') axesC.findprop('YLim')],...
        'PropertyPostSet',@localLimitListener);      
    end
    
    hFig = get(imAxes,'Parent');
    setappdata(hFig,'KenAxeResizeListeners',LimListener);
  end
  
  
  function localLimitListener(~,event)
    % when called, rescale all fonts in image
    imAxes = event.AffectedObject;
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      ud = getappdata(imAxes, 'TextZoomable_UserData');
    else
      ud = get(imAxes,'UserData');
    end
    
    fs = getBestFontSize(imAxes);
    for ii = 1:length(ud.handles)
      set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
    end
  end
%pc



























% 
% 
% 
%     % create the text
%     txtHandle = text(x,y,varargin{:});
% 
%     % detect its size relative to the figure, and set up listeners to resize
%     % it as the figure resizes, or axis limits are changed.
%     hAx = gca;
%     hFig = get(hAx,'Parent');
%     
%     fs = get(txtHandle, 'FontSize');
%     ratios = fs * diff(get(hAx,'YLim')) / max(get(hFig,'Position') .* [0 0 0 1]);
%     
%     
%     % append the handles and ratios to the user data - repeated calls will
%     % add each block of text to the list
%     ud = get(hAx, 'UserData');
%     if isfield(ud, 'ratios')
%         ud.ratios = [ud.ratios(:); ratios];
%     else
%         ud.ratios = ratios;
%     end
%     if isfield(ud, 'handles')
%         ud.handles = [ud.handles(:); txtHandle];
%     else
%         ud.handles = txtHandle;
%     end
%     
%     set(hAx,'UserData', ud);
%     localSetupPositionListener(hFig,hAx);
%     localSetupLimitListener(hAx);
% 
% end
% 
% 
% %% Helper Functions
% 
% function fs = getBestFontSize(imAxes)
%     % Try to keep font size reasonable for text
%     hFig = get(imAxes,'Parent');
%     hFigFactor = max(get(hFig,'Position') .* [0 0 0 1]);  
%     axHeight = diff(get(imAxes,'YLim'));
%     ud = get(imAxes,'UserData');  % stored in teh first user data.
%     fs = round(ud.ratios * hFigFactor / axHeight);    
%     fs = max(fs, 3);
% end
% 
% function localSetupPositionListener(hFig,imAxes)
%     % helper function to sets up listeners for resizing, so we can detect if
%     % we would need to change the fontsize
%     PostPositionListener = handle.listener(hFig,'ResizeEvent',...
%         {@localPostPositionListener,imAxes});
%     setappdata(hFig,'KenFigResizeListeners',PostPositionListener);
% end
% 
% function localPostPositionListener(~,~,imAxes) 
%     % when called, rescale all fonts in image
%     ud = get(imAxes,'UserData');
%     fs = getBestFontSize(imAxes);
%     for ii = 1:length(ud.handles)
%         set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
%     end   
% end
% 
% function localSetupLimitListener(imAxes)
%     % helper function to sets up listeners for zooming, so we can detect if
%     % we would need to change the fontsize
%     hgp     = findpackage('hg');
%     axesC   = findclass(hgp,'axes');
%     LimListener = handle.listener(imAxes,[axesC.findprop('XLim') axesC.findprop('YLim')],...
%         'PropertyPostSet',@localLimitListener);
%     hFig = get(imAxes,'Parent');
%     setappdata(hFig,'KenAxeResizeListeners',LimListener);
% end
% 
% function localLimitListener(~,event)
%     % when called, rescale all fonts in image
%     imAxes = event.AffectedObject;
%     ud = get(imAxes,'UserData');
%     fs = getBestFontSize(imAxes);
%     for ii = 1:length(ud.handles)
%         set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
%     end
% end
