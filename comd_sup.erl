-module(comd_sup).
-behaviour(supervisor).

%% API
-export([start_link/2]).

%% Supervisor callbacks
-export([init/1]).

-define(SERVER, ?MODULE).
-define(CMD, "./comd -x 123.0 -y 234.0 -z 345.0").


%%====================================================================
%% API
%%====================================================================

start_link(Pid, N) ->
    supervisor:start_link(?MODULE, [Pid, N]).
    %% supervisor:start_link({local, ?SERVER}, ?MODULE, [Pid, N]).

%%====================================================================
%% Supervisor callbacks
%%====================================================================

init([Pid, N]) ->
    L = lists:seq(0, N-1, 1),
    Children = lists:map(
        fun(I) ->
            {mk_name("tag_", I),
             {comd_srv, start_link, [mk_name("comd_", I), ?CMD, I, 0, Pid]},
             transient, 2000, worker, dynamic}
        end,
        L
    ),
    {ok, { {one_for_one, 3, 10}, Children }}.

%%====================================================================
%% Internal functions
%%====================================================================

mk_name(Base, N) ->
    list_to_atom(Base ++ integer_to_list(N)).



