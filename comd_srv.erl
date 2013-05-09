-module(comd_srv).
-behaviour(gen_server).

%% Server API
-export([start_link/5, stop/1]).

%% Client API
-export([set_stress/2]).

%% gen_server callbacks
-export([init/1, handle_call/3, handle_cast/2,
         handle_info/2, terminate/2, code_change/3]).

-define(SERVER(Cell), list_to_atom("comd_" ++ integer_to_list(Cell))).

-record(state, {name, command_line, index, time_step, port, hmm_pid}).

%%====================================================================
%% Server API
%%====================================================================

start_link(Name, CommandLine, Index, TimeStep, Pid) ->
    Result = gen_server:start_link(
               {local, ?SERVER(Index)},
               ?MODULE, [Name, CommandLine, Index, TimeStep, Pid], []),
    %io:format("~p: started~n", [Name]),
    Result.

stop(Name) ->
    gen_server:cast(?SERVER(Name), shutdown).

%%====================================================================
%% Client API
%%====================================================================

set_stress(Name, Stress) ->
    gen_server:cast(Name, {set_stress, Stress}).

%%====================================================================
%% gen_server callbacks
%%====================================================================

init([Name, CommandLine, Index, TimeStep, Pid]) ->
    io:format("+"),
    process_flag(trap_exit, true),
    Port = open_port({spawn, CommandLine}, [use_stdio, exit_status]),
    % Send our name and cell to HMM so it can tell us what to do.
    Pid ! {comd_started, {Index, Name}},
    {ok, #state{name = Name,
                command_line = CommandLine,
                index = Index,
                time_step = TimeStep,
                port = Port,
                hmm_pid=Pid}}.

handle_cast({set_stress, Stress},
            State = #state{name = Name, time_step = TimeStep, port = Port}) ->
    io:format("Server[~p] got stress: ~p~n", [Name, Stress]),
    [StressX, StressY] = Stress,
    Formatted = io_lib:format("~.6f ~.6f~n", [StressX, StressY]),
    Payload = list_to_binary(lists:flatten(Formatted)),
    io:format("[~p] ~p~n", [Port, Payload]),
    erlang:port_command(Port, Payload),
    {noreply, State#state{time_step = TimeStep + 1}}.

handle_call(_Message, _From, State) ->
    {reply, ok, State}.

%% Get replies from the external code. Messages preceded
%% by code that specifies function to be carried out here.
%%
%%   stress: resultant strain from comd
handle_info({Port, {data, Data}}, 
            State = #state{name = Name,
                           command_line = CommandLine,
                           index = Index,
                           time_step = TimeStep,
                           port = Port,
                           hmm_pid=Pid}) ->
    %%L = string:tokens(Data, " "),
    %%[CodeString|_] = L,
    %%Code = list_to_atom(CodeString),
    %Fun = fun string_to_num/1,
    %lists:foreach(Fun, L),
    %io:format("[~p]{~p} Data: ~p~n", [Name, TimeStep, Code]),
    %%process_code(Code, L, Index, TimeStep, Pid),
    Pid ! {Data},
    {noreply, State};
handle_info({Port, {exit_status, 234}}, State) ->
    % io:format("NORMAL: ~p ~p~n", [Port, 234]),
    {stop, normal, State};
handle_info({Port, {exit_status, 143}}, State) ->
    % io:format("KILLED: ~p ~p~n", [Port, 143]),
    exit(crashed_or_killed);
handle_info({'EXIT', Port, Status}, State) ->
    {noreply, State}.

terminate(normal, #state{name = Name}) ->
    %io:format("~p: shutdown cleanly (voluntary)~n", [Name]),
    ok;
terminate(shutdown, #state{name = Name}) ->
    io:format("~p: shutdown cleanly (forced)~n", [Name]),
    ok;
terminate(Reason, #state{name = Name}) ->
    io:format("~p: terminated because ~1024p~n", [Name, Reason]),
    ok.

code_change(_OldVsn, State, _Extra) ->
    {ok, State}.

%%====================================================================
%% Internal functions
%%====================================================================

string_to_num(S) ->
    case string:to_float(S) of
    {error,no_float} -> 
        list_to_integer(S);
     {F,_Rest} -> 
        % io:format("~p~n", [F]),
        F
    end.

%% These return message to be sent to hmm.
process_code(Code, Tokens, Ix, TS, Pid) when Code == strain ->
    io:format("Got strain: ~p~n", [Tokens]),
    TokensNoCode = lists:delete("strain", Tokens),
    Strain = lists:map(fun(E) -> string_to_num(E) end, TokensNoCode),
    Pid ! {strain, Ix, Strain}.
