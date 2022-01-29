%% clear vars 
sca;
close all;
clearvars;
tic; 
AssertOpenGL; 
%% experiment set-up
EXP.exp_vers = 'SCC_v1_part1'; 
EXP.matlab_vers = version; 
EXP.PTB_vers = PsychtoolboxVersion; 
if ispc
    EXP.os = [system_dependent('getos'), '_', system_dependent('getwinsys')]; 
else
    EXP.os = [system_dependent('getos')]; 
    Screen('Preference','SkipSyncTests', 1); 
end
EXP.start_time = datetime ('now');
EXP.keypad = false;
EXP.sub_num = input('Input subject number and press enter: ');
if isempty(dir(['data/*part1*' num2str(EXP.sub_num) '.mat']))
   EXP.save_name = ['data/' EXP.exp_vers '_sub_' num2str(EXP.sub_num) '.mat']; 
else 
   error('A data file exists for this subject number!'); 
end     
%% stimulus and screen parameters
% stimulus parameters
scene_size.x = 1000; 
scene_size.y = 800; 
target_size = 40; 
target_bold = target_size/7; 
target_color = [196 196 196]; 

% get screen measures 
screenNumber = max(Screen('Screens'));
screenx.w = WhiteIndex(screenNumber);
screenx.b = BlackIndex(screenNumber);
screenx.g = screenx.w / 2;
[w, wRect]=Screen('OpenWindow', screenNumber, screenx.b);
[screen_size.x, screen_size.y] = RectSize(wRect);
[screen_center.x, screen_center.y] = RectCenter(wRect);
% set priority for script execution to realtime priority
priorityLevel=MaxPriority(w);
Priority(priorityLevel);

% scene position
scene_pos = [screen_center.x - scene_size.x/2, ...
             screen_center.y - scene_size.y/2, ...
             screen_center.x + scene_size.x/2, ...
             screen_center.y + scene_size.y/2];  

% prepare fixation cross
fix_size = 50; 
fixCr = zeros(fix_size,fix_size);
fixCr(fix_size/2-2: fix_size/2+2, :) = unique(target_color);
fixCr(:, fix_size/2-2: fix_size/2+2) =  unique(target_color);  
fixposition = [screen_center.x-fix_size/2, screen_center.y-fix_size/2,...
               screen_center.x+fix_size/2, screen_center.y+fix_size/2];

%% timing parameters
dur_initblank = 1000/1000; % blank screen after the instruction 
dur_fix = 500/1000; 
dur_ITI = 750/1000; 
dur_rtlimit = 15000/1000; 
dur_break = 60; 
dur_beep = .01; 
%% sound params
beep_corr = 400; 
beep_incorr = 3000; 
%% Configure keyboard for response input
KbName('UnifyKeyNames');
if EXP.keypad
    key_left = KbName('1');        
    key_right = KbName('2');
else
    key_left = KbName('LeftArrow');        
    key_right = KbName('RightArrow');
    key_space = KbName('space'); 
end
%% Prepare instruction texts
ShowText.Instruction ='Find a letter T rotated 90 degrees in the scene.\n\nUsing arrow keys, indicate which direction the bottom of the T is pointing.\n\nRespond as quickly as possible.\n\nFind experimenter if you have questions.\n\nPress the space bar to proceed.'; 
ShowText.Break ='Take a break.\n\n The experiment will resume in 1 minute.'; 
ShowText.Resume ='Press the space bar\n\nwhen you are ready to start the next block.'; 
ShowText.End ='You finished this part of the study.\n\nThank you for your time!'; 
%% Define the grid space and target locations
% Define grid space
Grid.ny = 8; 
Grid.nx = 10; 
Grid.dim = 100; 
Grid.cellidx = 1:Grid.ny*Grid.nx; 
Grid.yidx = floor((Grid.cellidx - 1)/Grid.nx)+1; 
Grid.xidx = mod((Grid.cellidx - 1), 10) + 1; 
Grid.ycoord = Grid.dim .* Grid.yidx + scene_pos(2) - Grid.dim/2; % center coordinate of each grid cell 
Grid.xcoord = Grid.dim .* Grid.xidx + scene_pos(1) - Grid.dim/2;
Grid.mask = ones(1, Grid.ny * Grid.nx); 
Grid.mask(Grid.yidx == 1 | Grid.yidx == Grid.ny) = 0; 
Grid.mask(Grid.xidx == 1 | Grid.xidx == Grid.nx) = 0; 
Grid.mask((Grid.yidx == Grid.ny/2 | Grid.yidx == Grid.ny/2+1) & (Grid.xidx == Grid.nx/2 | Grid.xidx == Grid.nx/2+1)) = 0; 
Grid.quadrant = zeros(size(Grid.cellidx)); 
Grid.quadrant(Grid.xidx <= Grid.nx/2 & Grid.yidx <= Grid.ny/2) = 1; 
Grid.quadrant(Grid.xidx > Grid.nx/2 & Grid.yidx <= Grid.ny/2) = 2; 
Grid.quadrant(Grid.xidx <= Grid.nx/2 & Grid.yidx > Grid.ny/2) = 3; 
Grid.quadrant(Grid.xidx > Grid.nx/2 & Grid.yidx > Grid.ny/2) = 4; 
% [Grid.cellidx; Grid.xidx; Grid.yidx;  Grid.quadrant]

% Define target locs
targ.seed_q1 = randsample (Grid.cellidx(find(Grid.quadrant == 1 & Grid.mask == 1)), 4); 
targ.seed_q2 = randsample (Grid.cellidx(find(Grid.quadrant == 2 & Grid.mask == 1)), 4); 
targ.p = [targ.seed_q1(1:2) targ.seed_q2(1:2)];
targ.np = [targ.seed_q1(3:4) targ.seed_q2(3:4)];
for i = 1: 1: 4 % get the mirrored coords from the other condition (to control the distance from the center) 
    targ.np = [targ.np find(Grid.xidx == (Grid.nx+1 - Grid.xidx(targ.p(i))) & Grid.yidx == (Grid.ny+1 - Grid.yidx(targ.p(i))))]; 
    targ.p = [targ.p find(Grid.xidx == (Grid.nx+1 - Grid.xidx(targ.np(i))) & Grid.yidx == (Grid.ny+1 - Grid.yidx(targ.np(i))))];  
end 
targ.p = Shuffle(targ.p); 
targ.np = Shuffle(targ.np); 
targ.quad_p = Grid.quadrant(Grid.cellidx(targ.p)); 
targ.quad_np = Grid.quadrant(Grid.cellidx(targ.np)); 
%% Define trial set  
% This is the default set with trial index unshuffled 
trial.nBlock = 16;
trial.nTrial = 16; 
trial.idx = [1:trial.nTrial]; 
trial.cond = floor((trial.idx-1) / (length(trial.idx)/2))+1; % 1 = p, 2 = np
trial.targdir = mod(floor((trial.idx-1) / (length(trial.idx)/4)), 2)+1; % 1 = L, 2 = R
trial.targquad = mod((trial.idx-1), (length(trial.idx)/4))+1; %
for i = 1: 1: trial.nTrial
    if trial.cond(i) == 1 % p condition 
       tmp = find(targ.quad_p == trial.targquad(i)); 
        if trial.targdir(i) == 1 % L target
           trial.targloc(i) = targ.p(tmp(1)); 
        else % R target 
           trial.targloc(i) = targ.p(tmp(2));
        end 
    else % np condition
       tmp = find(targ.quad_np == trial.targquad(i)); 
        if trial.targdir(i) == 1 % L target
           trial.targloc(i) = targ.np(tmp(1)); 
        else % R target 
           trial.targloc(i) = targ.np(tmp(2));
        end 
    end 
end 
scheme.p = sort(randsample(trial.idx, trial.nTrial/2), 'ascend'); 
scheme.np = setdiff(trial.idx, scheme.p); 
% targ.p 1 2 ... 8 --> pair with the same scene scheme 
% targ.np 1 2 ... 8 --> pair with different scene scheme each time 
trial.scheme(trial.cond == 1) = scheme.p; % pre-determined
trial.scheme(trial.cond == 2) = 99; % scheme.np shuffled each block.    

% [trial.idx' trial.cond' trial.targdir' trial.targquad' trial.targloc' trial.scheme']
%% initialize datalog 
datalog.sceneorder = Shuffle([1:trial.nBlock]); % use this instead of iBlock so that the order of the scenes presented is not the same in all subjs 
datalog.schemecon = zeros(trial.nBlock, trial.nTrial); 
datalog.schemeidx = zeros(trial.nBlock, trial.nTrial); 
datalog.trial_idx = zeros(trial.nBlock, trial.nTrial); 
datalog.targquad = zeros(trial.nBlock, trial.nTrial); 
datalog.targloc = zeros(trial.nBlock, trial.nTrial); 
datalog.targdir = zeros(trial.nBlock, trial.nTrial); 
datalog.trialonset_stamp = zeros(trial.nBlock, trial.nTrial); 
datalog.stimon_stamp = zeros(trial.nBlock, trial.nTrial); 
datalog.resp_stamp = zeros(trial.nBlock, trial.nTrial); 
datalog.rt = zeros(trial.nBlock, trial.nTrial); 
datalog.acc = zeros(trial.nBlock, trial.nTrial); 
datalog.keypress = zeros(trial.nBlock, trial.nTrial); 
%% START THE EXPERIMENT 
RestrictKeysForKbCheck([key_left,key_right, key_space]);
HideCursor; 
targdirs = {'L', 'R'}; % correspond to 1 & 2

for iBlock = 1: 1: trial.nBlock 
    % get the randomized scene idx
    sceneidx = datalog.sceneorder(iBlock); 
    
    % randomize np schemes, targdir, and  trialidx within blocks 
    trial.scheme(trial.cond == 2) = Shuffle(scheme.np); % scheme.np shuffled each block. 
    trial_randidx = Shuffle(trial.idx); 
    trial.targdir(trial.cond == 1) = Shuffle(trial.targdir(trial.cond == 1)); 
    trial.targdir(trial.cond == 2) = Shuffle(trial.targdir(trial.cond == 2)); 
    
    % generate target stims for the trial set 
    for iTrial = 1: 1: trial.nTrial
        targ_pos(iTrial, :)  = [Grid.xcoord(trial.targloc(iTrial))-target_size/2, Grid.ycoord(trial.targloc(iTrial))-target_size/2, ...
                    Grid.xcoord(trial.targloc(iTrial))+target_size/2, Grid.ycoord(trial.targloc(iTrial))+target_size/2]; 
        targrectH(iTrial, :)  = [targ_pos(iTrial, 1), (targ_pos(iTrial, 2)+targ_pos(iTrial, 4))/2 - target_bold/2, ... 
            targ_pos(iTrial, 3) (targ_pos(iTrial, 2)+targ_pos(iTrial, 4))/2 + target_bold/2];  
        if targdirs{trial.targdir(iTrial)} == 'L'    
            targrectV(iTrial, :)  = [targ_pos(iTrial, 3)-target_bold, targ_pos(iTrial, 2), ...
                                     targ_pos(iTrial, 3) targ_pos(iTrial, 4)];  
        elseif targdirs{trial.targdir(iTrial)} == 'R'
            targrectV(iTrial, :)  = [targ_pos(iTrial, 1: 2), ...
                                     targ_pos(iTrial, 1)+target_bold targ_pos(iTrial, 4)];  
        else 
            warning('check the target direction'); 
        end 
    end    
    
    % show instruction message
    if mod(iBlock, 4) == 1 
        Screen('TextSize', w, 30);
        DrawFormattedText(w, ShowText.Instruction, 'center','center', WhiteIndex(w)); % prepare message
        Screen('Flip', w); % update the display to show the study message
        KbWait; % wait for kb input:

        % update the display to show blank screen for 1.5 seconds
        Screen('Flip', w);
        WaitSecs(dur_initblank);    
    end 
    
    for iTrial = 1: 1: trial.nTrial
        % replace iTrial with tidx = trial_randidx(iTrial) once you're ready!! 
        tidx = trial_randidx(iTrial); % tidx = iTrial; 
        
        % draw fixation cross
        fixcross = Screen('MakeTexture', w, fixCr);
        Screen('DrawTexture', w, fixcross, [], fixposition);
        [VBLTimestamp, FixOnsetTime] = Screen('Flip', w);
        WaitSecs(dur_fix);       
                 
        % stimulus presentation 
        scene_stim = imread(['v1_part1_stims/scheme_', num2str(trial.scheme(tidx)), '/', num2str(sceneidx), '.jpg']);
        scene_img = Screen('MakeTexture', w, scene_stim);
        Screen('DrawTextures', w, [scene_img], [], [scene_pos']);
        Screen('FillRect', w, target_color', [targrectH(tidx, :); targrectV(tidx, :)]'); 
        [VBLTimestamp, StimulusOnsetTime] = Screen('Flip', w);
 
        % get response 
        KbReleaseWait; % waits until all keys on the keyboard are released.
        KeyCode = zeros(1, 256); 
        KeyIsDown = 0; 
        ResponseTime = 0;    
        timeout = false; 
        
        while ~timeout
            [KeyIsDown, TimeStamp, KeyCode]=KbCheck;
            if KeyIsDown 
                if KeyCode(key_left) && (targdirs{trial.targdir(tidx)} == 'L')
                    datalog.acc(iBlock, iTrial) = 1; 
                    Beeper(beep_corr, dur_beep);
                elseif KeyCode(key_right) && (targdirs{trial.targdir(tidx)} == 'R')
                    datalog.acc(iBlock, iTrial) = 1; 
                    Beeper(beep_corr, dur_beep);
                else 
                    datalog.acc(iBlock, iTrial) = 0; 
                    Beeper(beep_incorr, dur_beep);
                end                
                break; 
            end 
            if GetSecs-StimulusOnsetTime > dur_rtlimit
               timeout = true; Beeper(beep_incorr, dur_beep);
            end 
        end 
        
        datalog.trial_idx(iBlock, iTrial) = tidx; 
        datalog.rt(iBlock, iTrial) = TimeStamp - StimulusOnsetTime;         
        datalog.schemecon(iBlock, iTrial) = trial.cond(tidx); 
        datalog.schemeidx(iBlock, iTrial) = trial.scheme(tidx); 
        datalog.sceneidx(iBlock, iTrial) = sceneidx; 
        datalog.targquad(iBlock, iTrial) = trial.targquad(tidx); 
        datalog.targloc(iBlock, iTrial) = trial.targloc(tidx); 
        datalog.targdir(iBlock, iTrial) = trial.targdir(tidx); 
        datalog.trialonset_stamp(iBlock, iTrial) = FixOnsetTime; 
        datalog.stimon_stamp(iBlock, iTrial) = StimulusOnsetTime; 
        datalog.resp_stamp(iBlock, iTrial) = TimeStamp;   
        if KeyCode(key_left)
            datalog.keypress(iBlock, iTrial) = 1;   
        elseif KeyCode(key_right)
            datalog.keypress(iBlock, iTrial) = 2;   
        end             
        
        % blank screen iti 
        Screen('Flip', w); 
        WaitSecs(dur_ITI);              
    end 
       
    if mod(iBlock, 4) == 0
        if iBlock ~= trial.nBlock 
            % show break message
            Screen('TextSize', w, 30);
            DrawFormattedText(w, ShowText.Break, 'center','center', WhiteIndex(w)); 
            Screen('Flip', w); 
            WaitSecs(dur_break);  

            % show resume message
            DrawFormattedText(w, ShowText.Resume, 'center','center', WhiteIndex(w));
            Screen('Flip', w);
            KbWait; % wait for kb input:
            
            KbReleaseWait;
            
        else % show exit message at end of the last block  
            Screen('TextSize', w, 30);
            DrawFormattedText(w, ShowText.End, 'center','center', WhiteIndex(w)); 
            Screen('Flip', w); 
            KbWait; sca; 
        end 
    end 
end

save (EXP.save_name, 'EXP', 'trial', 'scheme', 'datalog'); 

%% cleanup at end of experiment; switch Matlab/Octave back to priority 0 -- normal priority:
Screen('CloseAll');
ShowCursor;
fclose('all');
Priority(0);
toc;




