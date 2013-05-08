-module(hmm).

-export([run/1]).

run(N) ->
    {ok, PidCoMD} = comd_sup:start_link(self(), N),
    {ok, PidHMM}  = hmm_sup:start_link(),
    hmm_srv:send_stress(losolve, [23.9, 22.3]),
    unlink(PidCoMD),   % so shell doesn't crash
    unlink(PidHMM),    % so shell doesn't crash
    Servers = get_servers([], N),
    io:format("~p servers started: ~p~n", [length(Servers), Servers]),
    timesteps(1, Servers).

%% Returns the list of stated CoMD servers.  This list is
%% collected by accepting startup messages from each individual
%% server.
get_servers(L, N) when N > 0 ->
    receive
        {comd_started, {Cell, Name}} ->
            ok
    end,
    get_servers([{Cell, Name} | L], N-1);
get_servers(L, 0) ->
    lists:reverse(L).
    
%% Run the servers through a number of timesteps.
%% This is the main execution loop of the code.
%%
%%  1.  Start servers running by sending them strain.  This is the HO task.
%%  2.  Loop waiting for all servers to return their stress.
%%  3.  Calculate new strains based on the returned stresses. This is the LO task.
%%  4.  Repeat.
timesteps(N, S) when N > 0 ->
    lists:foreach(fun send_stress/1, S),
    Strain = get_strain([], length(S)),
    io:format("strains: ~p~n", [lists:sort(Strain)]),
    %update_stress(AllNames),
    timesteps(N-1, S);
timesteps(0, S) ->
    ok.

send_stress(E) ->
    {_, Server} = E,
    comd_srv:set_stress(Server, [-9863.09, 6353.2]).

get_strain(L, N) when N > 0 ->
    receive
        {strain, Cell, Strain} ->
            ok
    end,
    get_strain([{Cell, Strain} | L], N-1);
get_strain(L, 0) ->
    L.
