-module(hmm_sup).
-behaviour(supervisor).

%% API
-export([start_link/0]).

%% Supervisor callbacks
-export([init/1]).

-define(SERVER, ?MODULE).
-define(CMD, "./hmm").


%%====================================================================
%% API
%%====================================================================

start_link() ->
    supervisor:start_link(?MODULE, []).
    %% supervisor:start_link({local, ?SERVER}, ?MODULE, []).

%%====================================================================
%% Supervisor callbacks
%%====================================================================

init([]) ->
    HMM =
        {hmm,
         {hmm_srv, start_link, [losolve, ?CMD]},
         transient, 2000, worker, dynamic},

    Children = [HMM],
    {ok, { {one_for_one, 3, 10}, Children }}.

%%====================================================================
%% Internal functions
%%====================================================================
