-module(hmm_srv).
-behaviour(gen_server).

%% Server API
-export([start_link/2, stop/1]).

%% Client API
-export([echo/1, send_stress/2]).

%% gen_server callbacks
-export([init/1, handle_call/3, handle_cast/2,
         handle_info/2, terminate/2, code_change/3]).

-record(state, {name, command_line, port}).

%%====================================================================
%% Server API
%%====================================================================

start_link(Name, CommandLine) ->
    Result = gen_server:start_link(
               {local, Name},
               ?MODULE, [Name, CommandLine], []),
    %io:format("~p: started~n", [Name]),
    io:format("Result: ~p~n", [Result]),
    Result.

stop(Name) ->
    gen_server:cast(Name, shutdown).

%%====================================================================
%% Client API
%%====================================================================

echo(Name) ->
    gen_server:cast(Name, echo_message).

send_stress(Name, Stress) ->
    gen_server:call(Name, {send_stress, Stress}).

%%====================================================================
%% gen_server callbacks
%%====================================================================

init([Name, CommandLine]) ->
    io:format("Name: ~p   CMD: ~p~n", [Name, CommandLine]),
    process_flag(trap_exit, true),
    Port = open_port({spawn, CommandLine}, [use_stdio, exit_status]),
    % Send our name and cell to HMM so it can tell us what to do.
    {ok, #state{name = Name,
                command_line = CommandLine, port = Port}}.

handle_cast(echo_message,
            State = #state{name = Name,
                           command_line = CommandLine}) ->
    io:format("~p: ~p~n", [Name, CommandLine]),
    {noreply, State}.
%handle_cast({send_stress, Stress},
%            State = #state{name = Name, port = Port}) ->
%    [StressX, StressY] = Stress,
%    Formatted = io_lib:format("~.6f ~.6f~n", [StressX, StressY]),
%    Payload = list_to_binary(lists:flatten(Formatted)),
%    erlang:port_command(Port, Payload),
%    {noreply, State}.

handle_call({send_stress, Stress}, _From,
            State = #state{name = Name, port = Port}) ->
    [StressList| _] = Stress,
    send_to_port(Port, StressList),
    {reply, ok, State};
handle_call(_Message, _From, State) ->
    {reply, ok, State}.

%% Get replies from the external code. Messages preceded
%% by code that specifies function to be carried out here.
handle_info({Port, {exit_status, 234}}, State) ->
    % io:format("NORMAL: ~p ~p~n", [Port, 234]),
    {stop, normal, State};
handle_info({Port, {exit_status, 143}}, State) ->
    % io:format("KILLED: ~p ~p~n", [Port, 143]),
    exit(crashed_or_killed);
handle_info({'EXIT', Port, Status}, State) ->
    {noreply, State}.

terminate(normal, #state{name = Name}) ->
    io:format("~p: shutdown cleanly (voluntary)~n", [Name]),
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

send_to_port(_P, []) ->
    ok;
send_to_port(P, L) ->
    [F | R] = L,
    %io:format("L: ~p~nF: ~p~nR: ~p~n", [L,F,R]),
    Formatted = io_lib:format("~.6f", [F]),
    Payload = list_to_binary(lists:flatten(Formatted)),
    erlang:port_command(P, Payload),
    send_to_port(P, R).





