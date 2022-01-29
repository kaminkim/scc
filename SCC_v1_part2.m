%% clear vars 
sca;
close all;
clearvars;
tic; 
AssertOpenGL; 
%% experiment set-up
EXP.exp_vers = 'SCC_v1_part2'; 
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
if isempty(dir(['data/*part2*' num2str(EXP.sub_num) '.mat']))
   EXP.save_name = ['data/' EXP.exp_vers '_sub_' num2str(EXP.sub_num) '.mat']; 
else 
   error('A data file exists for this subject number!'); 
end     
%% load scene images 
stim_dir = 'v1_part2_stims/'; 
dlist = dir(stim_dir); 
dlist = dlist(arrayfun(@(x) x.isdir(1), dlist) == 1);
dlist = dlist(arrayfun(@(x) x.name(1), dlist) ~= '.');

sceneidx = 1; 
scenes = cell(384, 1); 
for iDir = 1: length(dlist)
    current_dir = [stim_dir dlist(iDir).name '/']; 
    flist = dir(current_dir); 
    flist = flist(arrayfun(@(x) x.isdir(1), flist) == 0);
    flist = flist(arrayfun(@(x) x.name(1), flist) ~= '.');
    for iFile = 1: length(flist)
        fname = [current_dir flist(iFile).name];  
        scenes{sceneidx} = fname; 
        sceneidx = sceneidx + 1; 
    end 
    clear flist 
end 

%% stimulus and screen parameters
% stimulus parameters
scene_size.x = 1000*0.8; 
scene_size.y = 800*0.8; 

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

% scale position (y) 
scale_pos = screen_center.y + (scene_size.y/2)*1.05; 

%% timing parameters
dur_initblank = 1000/1000; % blank screen after the instruction 
dur_ITI = 750/1000; 
dur_break = 60; 
%% Configure keyboard for response input
KbName('UnifyKeyNames');
if EXP.keypad
%     key_left = KbName('1');        
%     key_right = KbName('2');
else
    key_sure_new = KbName('1!'); 
    key_prob_new = KbName('2@'); 
    key_guess = KbName('3#'); 
    key_prob_old = KbName('4$'); 
    key_sure_old = KbName('5%'); 
    key_space = KbName('space'); 
end
%% Prepare instruction texts
ShowText.Instruction ='Please rate how confident you are \n\nthat the scene was presented in the first part of the experiment.\n\nPress the space bar to proceed.'; 
ShowText.Break ='Take a break.\n\n The experiment will resume in 1 minute.'; 
ShowText.Resume ='Press the space bar\n\nwhen you are ready to start the next block.'; 
ShowText.End ='You finished this part of the study.\n\nThank you for your time!'; 
ShowText.scalenum =    '1 ------- 2 ------- 3 ------- 4 ------- 5'; 
ShowText.scaletxt ='sure new                                   sure old'; 

%% initialize datalog 
datalog.scene = cell(size(scenes, 1), 1); 
datalog.resp = zeros(size(scenes, 1), 1); 
datalog.rt = zeros(size(scenes, 1), 1); 

%% START THE EXPERIMENT 
RestrictKeysForKbCheck([key_sure_new, key_prob_new, key_guess, key_prob_old, key_sure_old, key_space]);
HideCursor; 
 
% get the randomized scene idx
sceneidx = Shuffle([1:size(scenes, 1)])'; 
blocksize = length(sceneidx)/4; 
for iTrial = 1: 1: size(scenes, 1)
    
    % show instruction message
    if mod(iTrial, blocksize) == 1 
        Screen('TextSize', w, 30);
        DrawFormattedText(w, ShowText.Instruction, 'center','center', WhiteIndex(w)); % prepare message
        Screen('Flip', w); % update the display to show the study message
        KbWait; % wait for kb input:

        % update the display to show blank screen for 1.5 seconds
        Screen('Flip', w);
        WaitSecs(dur_initblank);    
    end 
    
    % stimulus presentation 
    present_scene = scenes{sceneidx(iTrial)} ; 
    scene_stim = imread(present_scene);
    scene_img = Screen('MakeTexture', w, scene_stim);
    Screen('DrawTextures', w, [scene_img], [], [scene_pos']);
    DrawFormattedText(w, ShowText.scalenum, 'center', scale_pos, WhiteIndex(w)); 
    DrawFormattedText(w, ShowText.scaletxt, 'center', scale_pos+25, WhiteIndex(w)); 
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
            if KeyCode(key_sure_new)
                datalog.resp(iTrial, 1)= 1; 
            elseif KeyCode(key_prob_new)
                datalog.resp(iTrial, 1)= 2; 
            elseif KeyCode(key_guess)
                datalog.resp(iTrial, 1)= 3; 
            elseif KeyCode(key_prob_old)
                datalog.resp(iTrial, 1)= 4; 
            elseif KeyCode(key_sure_old)
                datalog.resp(iTrial, 1)= 5; 
            end              
            break; 
        end 
    end 
  
    datalog.scene{iTrial} = present_scene;
    datalog.rt(iTrial, 1)    = TimeStamp - StimulusOnsetTime;               

    % blank screen iti 
    Screen('Flip', w); 
    WaitSecs(dur_ITI);              
       
    if mod(iTrial, blocksize) == 0
        if iTrial ~= size(scenes, 1)
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

save (EXP.save_name, 'datalog'); 

%% cleanup at end of experiment; switch Matlab/Octave back to priority 0 -- normal priority:
Screen('CloseAll');
ShowCursor;
fclose('all');
Priority(0);
toc;
