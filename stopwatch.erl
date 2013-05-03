-module(stopwatch).
-behaviour(gen_server).

%% Server API
-export([start_link/0, stop/0]).

%% Client API
-export([start_timer/0, stop_timer/0, read_timer/0]).

%% gen_server callbacks
-export([init/1, handle_call/3, handle_cast/2, handle_info/2, terminate/2, code_change/3]).

-define(SERVER, stopwatch).

-record(state, {mode = stopped,
                start_time = undefined,
                stop_time = undefined}).

%%====================================================================
%% Server API
%%====================================================================

start_link() ->
    Result = gen_server:start_link({local, ?SERVER}, ?MODULE, [], []),
    io:format("Stopwatch server started~n"),
    Result.

stop() ->
    gen_server:cast(?SERVER, shutdown).

%%====================================================================
%% Client API
%%====================================================================

start_timer() ->
    gen_server:cast(?SERVER, start_timer).

stop_timer() ->
    gen_server:cast(?SERVER, stop_timer).

read_timer() ->
    gen_server:call(?SERVER, read_timer).

%%====================================================================
%% gen_server callbacks
%%====================================================================

init([]) ->
    {ok, #state{}}.

handle_cast(shutdown, State) ->
    {stop, normal, State};

handle_cast(start_timer, State) ->
    StartTime = now(),
    io:format("Stopwatch started~n"),
    {noreply, State#state{mode = running, start_time = StartTime}};

handle_cast(stop_timer, State) ->
    StopTime = now(),
    io:format("Stopwatch stopped~n"),
    {noreply, State#state{mode = stopped, stop_time = StopTime}}.

handle_call(read_timer, _From, State = #state{mode = running}) ->
    Now = now(),
    Diff = calc_diff(Now, State#state.start_time),
    {reply, Diff, State};

handle_call(read_timer, _From, State = #state{mode = stopped}) ->
    Diff = calc_diff(State#state.stop_time, State#state.start_time),
    {reply, Diff, State}.

handle_info(_Info, State) ->
    {noreply, State}.

terminate(_Reason, _State) ->
    io:format("Stopwatch server is shutting down~n"),
    ok.

code_change(_OldVsn, State, _Extra) ->
    {ok, State}.

%%====================================================================
%% Internal functions
%%====================================================================

calc_diff(T2, T1) ->
    DiffMicroSecs = timer:now_diff(T2, T1),
    round(DiffMicroSecs / 1000).