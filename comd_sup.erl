-module(comd_sup).
-behaviour(supervisor).

%% API
-export([start_link/0]).

%% Supervisor callbacks
-export([init/1]).

-define(SERVER, ?MODULE).
-define(COMMAND_LINE, "./comd -x 123.0 -y 234.0 -z 345.0").

%%====================================================================
%% API
%%====================================================================

start_link() ->
    supervisor:start_link({local, ?SERVER}, ?MODULE, []).

%%====================================================================
%% Supervisor callbacks
%%====================================================================

init([]) ->
    Butt =  % Starting one server fo now.
        {tag1,
         {comd_srv, start_link, [butt, ?COMMAND_LINE, 0, 0]},
         permanent, 2000, worker, dynamic},

    Children = [Butt],
    {ok, {{one_for_one, 3, 10}, Children}}.

%%====================================================================
%% Internal functions
%%====================================================================

