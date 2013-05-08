-module(hmm_sup).
-behaviour(supervisor).

%% API
-export([start_link/2]).

%% Supervisor callbacks
-export([init/1]).

-define(SERVER, ?MODULE).
-define(CMD, "./OldCoMD/LOsolve").


%%====================================================================
%% API
%%====================================================================

start_link(Pid, N) ->
    supervisor:start_link({local, ?SERVER}, ?MODULE, [Pid, N]).

%%====================================================================
%% Supervisor callbacks
%%====================================================================

init([Pid, N]) ->
    HMM =
        {hmm,
         {hmm_srv, start_link, []},
         transient, 2000, worker, dynamic},

    Children = [HMM],
    {ok, { {one_for_one, 3, 10}, Children }}.

%%====================================================================
%% Internal functions
%%====================================================================
