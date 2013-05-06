-module(comd_srv).
-behaviour(gen_server).

%% Server API
-export([start_link/4, stop/1]).

%% Client API
-export([echo/1]).

%% gen_server callbacks
-export([init/1, handle_call/3, handle_cast/2,
         handle_info/2, terminate/2, code_change/3]).

-define(SERVER(Name), list_to_atom("comd_srv_" ++ atom_to_list(Name))).

-record(state, {name, command_line, index, time_step}).

%%====================================================================
%% Server API
%%====================================================================

start_link(Name, CommandLine, Index, TimeStep) ->
    Result = gen_server:start_link(
               {local, ?SERVER(Name)},
               ?MODULE, [Name, CommandLine, Index, TimeStep], []),
    io:format("~p: started~n", [Name]),
    Result.

stop(Name) ->
    gen_server:cast(?SERVER(Name), shutdown).

%%====================================================================
%% Client API
%%====================================================================

echo(Name) ->
    gen_server:cast(?SERVER(Name), echo_message).

%%====================================================================
%% gen_server callbacks
%%====================================================================

init([Name, CommandLine, Index, TimeStep]) ->
    process_flag(trap_exit, true),
    Port = open_port({spawn, CommandLine}, [use_stdio, exit_status]),
    Payload = list_to_binary("1.0 2.0 3.0\n") ,
    % io:format("Opened the port: ~w~n", [Port]),
    erlang:port_command(Port, Payload),
    % io:format("Sent command to port: ~p~n", [Payload]),
    {ok, #state{name = Name,
                command_line = CommandLine,
                index = Index,
                time_step = TimeStep}}.

handle_cast(echo_message,
            State = #state{name = Name,
                           command_line = CommandLine,
                           index = Index,
                           time_step = TimeStep}) ->
    % io:format("~p: ~p ~p ~p~n", [Name, CommandLine, Index, TimeStep]),
    {noreply, State#state{time_step = TimeStep + 1}}.

handle_call(_Message, _From, State) ->
    {reply, ok, State}.

handle_info({Port, {data, Data}}, State) ->
    % io:format("Received data: ~w~n", [Data]),
    L = string:tokens(Data, " "),
    % io:format("L: ~p~n", [L]),
    Fun = fun string_to_num/1,
    lists:foreach(Fun, L),
    % io:format("Data: ~p~n", [Data]),
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

string_to_num(S) ->
    case string:to_float(S) of
    {error,no_float} -> 
        list_to_integer(S);
     {F,_Rest} -> 
        % io:format("~p~n", [F]),
        F
    end.