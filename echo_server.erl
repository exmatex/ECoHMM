-module(echo_server).
-behaviour(gen_server).

%% Server API
-export([start_link/3, stop/1]).

%% Client API
-export([echo/1, assassinate/1]).

%% gen_server callbacks
-export([init/1, handle_call/3, handle_cast/2,
         handle_info/2, terminate/2, code_change/3]).

-define(SERVER(Name), list_to_atom("echo_server_" ++ atom_to_list(Name))).

-record(state, {name, message, interval, count}).

%%====================================================================
%% Server API
%%====================================================================

start_link(Name, Message, Interval)
  when is_atom(Name), is_list(Message), is_integer(Interval) ->
    Result = gen_server:start_link(
               {local, ?SERVER(Name)},
               ?MODULE, [Name, Message, Interval], []),
    io:format("~p: started~n", [Name]),
    Result.

stop(Name) ->
    gen_server:cast(?SERVER(Name), shutdown).

%%====================================================================
%% Client API
%%====================================================================

echo(Name) ->
    gen_server:cast(?SERVER(Name), echo_message).

assassinate(Name) ->
    gen_server:cast(?SERVER(Name), die_horribly).

%%====================================================================
%% gen_server callbacks
%%====================================================================

init([Name, Message, Interval]) ->
    process_flag(trap_exit, true),
    timer:apply_interval(Interval, ?MODULE, echo, [Name]),
    {ok, #state{name = Name,
                message = Message,
                interval = Interval,
                count = 1}}.

handle_cast(shutdown, State) ->
    {stop, normal, State};

handle_cast(echo_message,
            State = #state{name = Name,
                           message = Message,
                           count = Count}) ->
    io:format("~p: ~s (#~B)~n", [Name, Message, Count]),
    {noreply, State#state{count = Count + 1}}.

handle_call(_Message, _From, State) ->
    {reply, ok, State}.

handle_info(_Info, State) ->
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