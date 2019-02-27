%Script to run PRISDIL the Prisoners' Dilemma simulation
t0=cputime;
player=[];
score=[];
tft=[];
opp_score=[];
%[player,score,opp_score]=prisdil(genes,pop_sz,n_o_mvs,...
%n_o_games,score,indiv_opp_score,mu,xovr)
[player,score,opp_score]=prisdil(70,20,150,20,score,...
opp_score,0.01,0.5);
disp('Total CPU time:');cputime-t0

function [player,score,opp_score]=prisdil(genes,pop_sz,n_o_mvs,n_o_games,scores,opp_score,mu,xovr)
%PRISDIL Driver of the prisoner's dilemma experiment
%Generate the population of random startegies or players.
%Each consisting of 5 columns of 70 genes.
%The player is therefore a matrix PLAYER[70,5]
%Column 1 contains the scores for the seventy possible outcoes
%Column 2 the 1st outcome
% " 3 2nd ""
% " 4 " 3rd ""
% " 5 contains the random (possible) moves for this
% player (the random strategy)
global tft_mvs adj_mvs evolv
%log_<file> will contain the results for each player and the opponent
log_score=[];
log_opp_score=[];
%pop will contain the population. It changes with every generation.
pop = [];
pop = genpop(pop_sz,genes);
%n_o_games games (or generations) are played
log_tft_mvs=[];
log_adj_mvs=[];
%Set when a mutation should take place; here every 2 generations;
mut_time=2;
evolv='true ';
t=cputime;
for game_n=1:n_o_games
game_n
%Operate a mutation if it is time
if game_n==mut_time & game_n ~=n_o_games
pop=mutate(pop,pop_sz,genes,mu);
mut_time=mut_time+2;
end;
%Each of the generated strategies or players plays n_o_games games
%of n_o_mvs moves with TIT_FOR_TAT and ADJUSTER
for indiv=1:pop_sz
player(1:genes,1:5)=pop(1:genes,(indiv-1)*5+1:(indiv-1)*5+5);
indiv_score=0;
indiv_opp_score=0;
%Player's opponent is TIT_FOR_TAT
tft_mvs=[];
[player,indiv_mvs,indiv_score,indiv_opp_score]=...
tit_for_tat(player,genes,n_o_mvs,indiv_mvs,indiv_score,...
indiv_opp_score);
% Store moves
log_tft_mvs=[log_tft_mvs;tft_mvs];
tft_mvs=[];
%Player's opponent is ADJUSTER
adj_mvs=[];
[player,indiv_mvs,indiv_score,indiv_opp_score]=...
adjuster(player,genes,n_o_mvs,indiv_mvs,indiv_score,...
indiv_opp_score);
% Store scores (Average over the two games)
score=[score indiv_score/2];
opp_score=[opp_score indiv_opp_score/2];
log_adj_mvs=[log_adj_mvs;adj_mvs];
adj_mvs=[];
% [player,indiv_score,tft_mvs,indiv_opp_score]=...
% evlv_play(indiv,player,genes,n_o_mvs,indiv_score,...
% tft_mvs,indiv_opp_score);
score=[];
opp_score=[];
end;
save pop;
disp('CPU time taken to evolve the strategy:'); cputime-t
evolv='false';
%Find evolved strategy with best average score over
%n_o_games and its rules
bst_player=[];
%[strat_score,bst_strat]=max(log_score(game_n,1:pop_sz))
[strat_score,bst_strat]=max(mean(log_score));
bst_player(1:genes,1:5)=pop(1:genes,...
(bst_strat-1)*5+1:(bst_strat-1)*5+5);
save bst_player;
%Output readably coded rules of strategy
coded_rul=[];
ruls_score=[];
[ruls_score,coded_rul]=...
code_rul(ruls_score,coded_rul,bst_player,genes);
save coded_rul;
%Play evolved strategy against TFT
disp('Do you want to see the performance of this evolved player against TFT?');
input('y or n between single quotes : ');
if ans=='y'
indiv_score=0;
indiv_opp_score=0;
[bst_player,indiv_mvs,indiv_score,indiv_opp_score]=...
tit_for_tat(bst_player,genes,n_o_mvs,indiv_mvs,indiv_score,...
indiv_opp_score);
disp('Player and TFT scores: ')
indiv_score
indiv_opp_score
% indiv_mvs
% tft_mvs
end
%Check if the strategy has been altered (It should not!)
coded_rul_aftr_tft=[];
ruls_score=[];
[ruls_score,coded_rul_aftr_tft]=...
code_rul(ruls_score,coded_rul_aftr_tft,bst_player,genes);
save coded_rul_aftr_tft;
%Play evolved strategy against ADJUSTER
disp('Do you want to see the performance of this evolved player against ADJUSTER?');
input('y or n between single quotes : ');
if ans=='y'
adj_mvs=[];
indiv_score=0;
indiv_opp_score=0;
[bst_player,indiv_mvs,indiv_score,indiv_opp_score]=...
adjuster(bst_player,genes,n_o_mvs,indiv_mvs,...
indiv_score,indiv_opp_score);
disp('Player and Adjuster scores: ')
indiv_score
indiv_opp_score
end
end
function init_ruls=start_ruls(genes)
%Builds the initial set of all possible outcomes
%from the num_mvs previous moves.
%Here 3 previous are only considered
%Second column
init_ruls=[];
for i=1:16
init_ruls(i,2)=3;
init_ruls(16+i,2)=5;
init_ruls(32+i,2)=0;
init_ruls(48+i,2)=1;
end;
for i=1:6
outcome=rand;
if rand>=.75
init_ruls(64+i,2)=3;
elseif outcome<.75 & outcome>=.5
init_ruls(64+i,2)=5;
elseif outcome<.5 & outcome>=2
init_ruls(64+i,2)=0;
else
init_ruls(64+i,2)=1;
end;
end;
%Third column
for i=1:4
init_ruls(i ,3)=3;
init_ruls(4+i ,3)=5;
init_ruls(8+i ,3)=0;
init_ruls(12+i,3)=1;
init_ruls(16+i,3)=3;
init_ruls(20+i,3)=5;
14
init_ruls(24+i,3)=0;
init_ruls(28+i,3)=1;
init_ruls(32+i,3)=3;
init_ruls(36+i,3)=5;
init_ruls(40+i,3)=0;
init_ruls(44+i,3)=1;
init_ruls(48+i,3)=3;
init_ruls(52+i,3)=5;
init_ruls(56+i,3)=0;
init_ruls(60+i,3)=1;
end;
for i=1:6
outcome=rand;
if rand>=.75
init_ruls(64+i,3)=3;
elseif outcome<.75 & outcome>=.5
init_ruls(64+i,3)=5;
elseif outcome<.5 & outcome>=2
init_ruls(64+i,3)=0;
else
init_ruls(64+i,3)=1;
end;
end;
%Fourth column
init_ruls(1,4)=3;
init_ruls(2,4)=5;
init_ruls(3,4)=0;
init_ruls(4,4)=1;
for j=2:16
init_ruls(4*(j-1)+1:4*(j-1)+4,4)=init_ruls(1:4,4);
end;
for i=1:6
outcome=rand;
if rand>=.75
init_ruls(64+i,3)=3;
elseif outcome<.75 & outcome>=.5
init_ruls(64+i,3)=5;
elseif outcome<.5 & outcome>=2
init_ruls(64+i,3)=0;
else
init_ruls(64+i,3)=1;
end;
end;
%First column contains the score of each rule randomly generrated
%ie:Before start we initialise player's rules scores to what they are
15
%initially, eg: if we have RRR for outcomes in rule 19 say then
%the score for that rule is initialised to 3*3, ie player(19,1)=3*3=9
for i=1:genes
init_ruls(i,1)=sum(init_ruls(i,2:4));
end;
function pop=genpop(pop_sz,genes)
%Generates a population of individuals (chromosomes) of length genes
%Columns 1,2,3,4 are the same for every player at start. A matrix
%'outcomes[genes,4]' is thus required to build up each player.
outcomes=start_ruls(genes);
for indiv=1:pop_sz
for num_gen=1:genes
if rand>0.5
pop(num_gen,(indiv-1)*5+5)=1;
else
pop(num_gen,(indiv-1)*5+5)=0;
end;
end;
pop(1:genes,(indiv-1)*5+1:(indiv-1)*5+4)=outcomes(1:genes,1:4);
end;
function pop=cross_ovr(pop,pop_sz,genes,log_score,game_n)
%Function CROSS_OVR operates a single point crossover on
%the best performing individuals so far
%Approach: find half pop_sz of individuals with best
%average score over games played so far; mate them to
%produce next generation. Keep the population size constant
%Implementation: sort individuals according to their
%average score against opponents;
%sorted scores are in srt_av_scr in ascending order and player
%numbers are in indiv_n.
if game_n==1
[srt_av_scr,indiv_n]=sort(log_score);
else
[srt_av_scr,indiv_n]=sort(mean(log_score(1:game_n,1:pop_sz)));
end
% Find the random pairs to mate
mat_pairs=[];
frst_prtnr=pop_sz;
half_pop=floor(pop_sz/2);
16
while frst_prtnr>=half_pop+1
scnd_prtnr=half_pop+ceil(half_pop*rand);
while frst_prtnr==scnd_prtnr
scnd_prtnr=half_pop+ceil(half_pop*rand);
end;
pair=[indiv_n(frst_prtnr) indiv_n(scnd_prtnr)];
mat_pairs=[mat_pairs;pair];
frst_prtnr=frst_prtnr-1;
end;
%mat_pairs
%Mating, single point cross_over
new_pop=[];
n_pairs=1;
indiv=1;
while n_pairs<=floor(pop_sz/2)
x_ovr_pnt=floor(genes*rand);
%First offspring including the whole description of player,
%ie 5 columns
new_pop(1:x_ovr_pnt,(indiv-1)*5+1:(indiv-1)*5+5)=...
pop(1:x_ovr_pnt,(mat_pairs(n_pairs,1)-1)*5+1:...
(mat_pairs(n_pairs,1)-1)*5+5);
new_pop(x_ovr_pnt+1:genes,(indiv-1)*5+1:(indiv-1)*5+5)=...
pop(x_ovr_pnt+1:genes,(mat_pairs(n_pairs,2)-1)*5+1:...
(mat_pairs(n_pairs,2)-1)*5+5);
%Second offspring
indiv=indiv+1;
new_pop(1:x_ovr_pnt,(indiv-1)*5+1:(indiv-1)*5+5)=...
pop(1:x_ovr_pnt,(mat_pairs(n_pairs,2)-1)*5+1:...
(mat_pairs(n_pairs,2)-1)*5+5);
new_pop(x_ovr_pnt+1:genes,(indiv-1)*5+1:(indiv-1)*5+5)=...
pop(x_ovr_pnt+1:genes,(mat_pairs(n_pairs,1)-1)*5+1:...
(mat_pairs(n_pairs,1)-1)*5+5);
indiv=indiv+1;
n_pairs=n_pairs+1;
end;
pop=new_pop;
end
function pop=mutate(pop,pop_sz,genes,mu)
%Function MUTATE that operates mutations on all chromosomes at
%regular intervals Let n_o_genes_t_mu be the number of genes to mutate
n_o_genes_t_mu=ceil(genes*mu);
for indiv=1:pop_sz
for count=1:n_o_genes_t_mu
gene_t_mu=ceil(genes*rand);
17
if pop(gene_t_mu,(indiv-1)*5+5) == 0
pop(gene_t_mu,(indiv-1)*5+5)=1;
else
pop(gene_t_mu,(indiv-1)*5+5)=0;
end
end
end
end
function [player,indiv_score,tft,indiv_opp_score]=evlv_play(indiv,player,genes,n_o_mvs,indiv_score,tft,indiv_opp_score)
%A game of n_o_mvs moves is played by individual indiv against tft.
%Their corresponding scores are recorded in array score(1:20)
%tft is a simple array of length n_o_mv tft[1:n_mv]
%corresponding to the n_mv moves of tft strategy
%tft plays first and it is a C, note that C is represented
%by 1 and D by 0 also R is represented by its value 3,
%T by 5, S by 0 and P by 1
indiv_score=0;
indiv_opp_score=0;
tft=[tft 1];
for mv_n=1:n_o_mvs
%Instead of choosing the best rule, we chose randomly
%among the 50% best ones
%Sort them in ascending order first, then chose
%randomly from the top half
[srt_ruls_scr,srt_ruls]=sort(player(1:genes,1));
half_genes=floor(genes/2);
bst_rul=srt_ruls(half_genes+ceil(half_genes*rand));
action=player(bst_rul,5);
% Update bst_rul of player
player(bst_rul,2)=player(bst_rul,3);
player(bst_rul,3)=player(bst_rul,4);
if tft(mv_n)==1 & action==1
% ie R
player(bst_rul,4)= 3;
indiv_opp_score=indiv_opp_score+3;
elseif tft(mv_n)==1 & action==0
% ie T
player(bst_rul,4)= 5;
indiv_opp_score=indiv_opp_score+0;
elseif tft(mv_n)==0 & action==1
% ie S
player(bst_rul,4)= 0;
indiv_opp_score=indiv_opp_score+5;
else
% ie P
player(bst_rul,4)= 1;
indiv_opp_score=indiv_opp_score+1;
end;
% Update score for that rule
player(bst_rul,1)=player(bst_rul,1)+player(bst_rul,4);
% Update score of player
indiv_score=indiv_score+player(bst_rul,4);
% update tft
tft=[tft action];
end;
end;
%Player with the TIT-FOR-TAT startegy
function [player,indiv_mvs,indiv_score,indiv_opp_score]=...
it_for_tat(player,genes,n_o_mvs,indiv_mvs,indiv_score,...
indiv_opp_score)
%A game of n_o_mvs moves is played by individual indiv
%against tft_mvs.
%Their corresponding scores are recorded in array score(1:20)
global tft_mvs evolv
%tft_mvs is a simple array of length n_o_mv tft_mvs[1:n_mv]
%corresponding to the n_mv moves of tft_mvs strategy
%tft_mvs plays first and it is a C, note that C is represented
%by 1 and D by 0 also R is represented by its value 3, T by 5,
%S by 0 and P by 1
%Set first the three preceding moves (more precisely outcomes) of TFT
%to random values
%These outcomes effectively make up a tft rule; call it tft_rul.
tft_rul=[];
for init_mv=1:3
outcm=rand;
if outcm<0.5,
if outcm<0.25,
tft_rul=[tft_rul 0];
else
tft_rul=[tft_rul 1];
end
elseif outcm<0.75,
19
tft_rul=[tft_rul 3] ;
else
tft_rul=[tft_rul 5];
end
end
n_o_t_bst_rul1_usd=0;
tft_mvs=[tft_mvs 1];
indiv_mvs=[];
score_so_far=[];
opp_score_so_far=[];
log_cntr=1;
for mv_n=1:n_o_mvs
%Instead of choosing among the 50% best rules, we look for the one
%which matches the last three moves of the game and apply it.
%This requires setting TFT to three random moves, as done in tft_rul.
%find rule matching that of tft, i.e what to do after a sequence of 3
%given moves eg: What to do after 333, corresponding to RRR? etc..
bst_rul=0;
for rule_n=1:genes
if player(rule_n,2:4)==tft_rul(1:3),
bst_rul=rule_n;
BST_RUL1=[bst_rul tft_rul(1:3)];
n_o_t_bst_rul1_usd=n_o_t_bst_rul1_usd+1;
break;
end
end
if bst_rul==0
%If no such sequence is available we look for what matches the
%last two moves
for rule_n=1:genes
if player(rule_n,3:4)==tft_rul(2:3),
bst_rul=rule_n;
BST_RUL2=[bst_rul tft_rul(2:3)];
break;
end
end
end
%If no such seq is available then look for a rule with the last
%outcome matching last outcome of tft_rul
if bst_rul==0
for rule_n=1:genes
if player(rule_n,4)==tft_rul(3),
bst_rul=rule_n;
20
BST_RUL3=[bst_rul tft_rul(3)];
break;
end
end
end
%Error trap: if bst_rul is still 0 then there is something
%wrong with player
if bst_rul==0
disp('bst_rul is still 0; It should not be. Check player');
bst_player
save player;
save tft_rul;
end;
%Otherwise
action=player(bst_rul,5);
% Update rule of tft (Note that this rule is updated in
% evolution satge and in performance play stage as well)
tft_rul(1)=tft_rul(2);
tft_rul(2)=tft_rul(3);
%
indiv_mvs=[indiv_mvs action];
if tft_mvs(mv_n)==1 & action==1
% ie outcome is R
mv_outcm=3;
tft_rul(3)=3;
elseif tft_mvs(mv_n)==1 & action==0
% ie outcome is T
mv_outcm=5;
tft_rul(3)=0;
elseif tft_mvs(mv_n)==0 & action==1
% ie outcome is S
mv_outcm=0;
tft_rul(3)=5;
else
% ie outcome is P
mv_outcm=1;
tft_rul(3)=1;
end
% Update bst_rul of player if in evolution stage
if evolv=='true ',
player(bst_rul,2)=player(bst_rul,3);
player(bst_rul,3)=player(bst_rul,4);
player(bst_rul,4)=mv_outcm;
end
% Update total score for the rule applied in this move
player(bst_rul,1)=player(bst_rul,1)+mv_outcm;
% Update total score of player and tft
indiv_score=indiv_score+mv_outcm;
indiv_opp_score=indiv_opp_score+tft_rul(3);
% update tft_mvs
tft_mvs=[tft_mvs action];
% Log scores every 10 moves for each player
log_cntr=log_cntr+1;
if log_cntr>10
score_so_far=[score_so_far indiv_score];
opp_score_so_far=[opp_score_so_far indiv_opp_score];
log_cntr=1;
end
end
if evolv=='true ', return; end;
tft_scores=[score_so_far' opp_score_so_far'];
save tft_scores;
plot(tft_scores);
title('Performance of evolved strategy against TFT');
xlabel('moves*10');
ylabel('scores');
end

%Player with strategy ADJUSTER
function [player,indiv_mvs,indiv_score,indiv_opp_score]=adjuster(player,genes,n_o_mvs,indiv_mvs,indiv_score,indiv_opp_score)
%A game of n_o_mvs moves is played by individual
%indiv against adj_mvs.
%Their corresponding scores are recorded in array score(1:20)
global adj_mvs evolv
%adj_mvs is a simple array of length n_o_mv adj_mvs[1:n_mv]
%corresponding to the n_mv moves of adj_mvs strategy
%ADJUSTER strategy relies on a high rate of defection, thus
%attempting to exploit the opponent. But, when the opponent
%is showing its teeth as well, ADJUSTER revises, i.e.
%'adjusts', its defection rate. In details adjuster plays as
%follows: Play 2D in every 3 successive moves. If oponent plays
%2D's in successive moves, then revise rate to just 1D in 3 moves.
%When opponent plays 2C's in a row then change rate to original
%value, ie 2D's in 3 successive moves.
%We set the 3 first outcomes to random ones corresponding to
%0 0 1. There are 8 possible sequences of 3 outcomes:
adj_rul=[];
rnd_strt_ruls=[1 1 3;1 5 3;1 5 0;1 1 1;5 1 3;5 1 0;5 5 0;5 5 3];
indx_o_rnd_rul=ceil(rand*8);
adj_rul=[adj_rul rnd_strt_ruls(indx_o_rnd_rul,1:3)];
%ADJUSTR plays first and it is a D, note that C is represented
%by 1 and D by 0
%also R is represented by its value 3, T by 5, S by 0 and P by 1
%We arbitrarily set the first cycle of moves to 0 0 1, i.e. DDC
adj_mvs=[adj_mvs 0 0 1];
indiv_mvs=[];
score_so_far=[];
opp_score_so_far=[];
log_cntr=1;
set_strat_rate_1=[1 0 1;1 1 0;0 1 1];%1 defection in evrey 3 moves
set_strat_rate_2=[0 0 1;0 1 0;1 0 0];%2 defections in evrey 3 moves
%At start rate of defection is 2
rate=2;
cycle=1;
change='true ';
for mv_n=1:n_o_mvs
%Instead of choosing among the 50% best rules, we look for the one
%which matches the last three moves of the game and apply it. This
%requires setting ADJ to three random moves, as done in adj_rul.
%find rule matching that of ADJ, i.e what to do when a sequence of 3
%given moves occurs. eg:What to do after moves 333,
%corresponding to RRR etc..
bst_rul=0;
for rule_n=1:genes
if player(rule_n,2:4)==adj_rul(1:3),
bst_rul=rule_n;
break;
end
end
if bst_rul==0
%If no such sequence is available we look for what matches
%the last two moves
for rule_n=1:genes
if player(rule_n,3:4)==adj_rul(2:3),
bst_rul=rule_n;
break;
end
end
end
23
%If no such seq is available then look for a rule with the
%last outcome matching last outcome of adj_rul
if bst_rul==0
for rule_n=1:genes
if player(rule_n,4)==adj_rul(3),
bst_rul=rule_n;
break;
end
end
end
%Error trap: if bst_rul is still 0 then there is something
%wrong with player
if bst_rul==0
disp('bst_rul is still 0; It should not be Check player'); bst_player
save player;
save adj_rul;
end;
%Otherwise
action=player(bst_rul,5);
indiv_mvs=[indiv_mvs action];
% update adj_mvs
if mv_n>1 & change=='true '
if action==0 & indiv_mvs(mv_n-1)==0,
rate=1;
change='false';
elseif action==1 & indiv_mvs(mv_n-1)==1,
rate=2;
change='false';
end
end
%
if cycle>3
if rate==2 & change=='false',
seq_o_mvs=ceil(3*rand);
adj_mvs=[adj_mvs set_strat_rate_2(seq_o_mvs,1:3)];
change='true ';
elseif rate==1 & change=='false'
seq_o_mvs=ceil(3*rand);
adj_mvs=[adj_mvs set_strat_rate_1(seq_o_mvs,1:3)];
change='true ';
else
seq_o_mvs=ceil(3*rand);
adj_mvs=[adj_mvs set_strat_rate_2(seq_o_mvs,1:3)];
rate=2;
24
end
cycle=1;
end
% Update rule of adj (Note that this rule is updated in
% evolution satge and in performance play stage as well)
adj_rul(1)=adj_rul(2);
adj_rul(2)=adj_rul(3);
%
if adj_mvs(mv_n)==1 & action==1
% ie R
mv_outcm= 3;
adj_rul(3)=3;
elseif adj_mvs(mv_n)==1 & action==0
% ie T
mv_outcm= 5;
adj_rul(3)=0;
elseif adj_mvs(mv_n)==0 & action==1
% ie S
mv_outcm= 0;
adj_rul(3)=5;
else
% ie P
mv_outcm= 1;
adj_rul(3)=1;
end
% Update bst_rul of player if in evolution stage
if evolv=='true ',
player(bst_rul,2)=player(bst_rul,3);
player(bst_rul,3)=player(bst_rul,4);
player(bst_rul,4)=mv_outcm;
end
% Update total score for the rule applied in this move
player(bst_rul,1)=player(bst_rul,1)+mv_outcm;
% Update total score of player
indiv_score=indiv_score+mv_outcm;
indiv_opp_score=indiv_opp_score+adj_rul(3);
cycle=cycle+1;
% Log scores every 10 moves for each player
log_cntr=log_cntr+1;
if log_cntr>10
score_so_far=[score_so_far indiv_score];
opp_score_so_far=[opp_score_so_far indiv_opp_score];
25
log_cntr=1;
end
end
if evolv=='true ', return; end;
score_so_far;
opp_score_so_far;
adj_scores=[score_so_far' opp_score_so_far'];
save adj_scores;
plot(adj_scores);
title('Performance of evolved strategy against ADJUSTER');
xlabel('moves*10');
ylabel('scores');
%indiv_mvs
%adj_mvs
end

%Other startegies can be include in the same way
function [ruls_score,coded_rul]=code_rul(ruls_score,coded_rul,bst_player,genes)
%Write strategy in a more readable form
%The scores of the rules are in array ruls_score
ruls_score(1:genes)=bst_player(1:genes,1);
for rul_n=1:genes
rule=[];
for outcome_n=2:4
if bst_player(rul_n,outcome_n)==3
rule=[rule 'R'];
elseif bst_player(rul_n,outcome_n)==5
rule=[rule 'T'];
elseif bst_player(rul_n,outcome_n)==0
rule=[rule 'S'];
else
rule=[rule 'P'];
end;
end;
if bst_player(rul_n,5)==0
rule=[rule ' D'];
else
rule=[rule ' C'];
end;
coded_rul=[coded_rul;rule];
end;