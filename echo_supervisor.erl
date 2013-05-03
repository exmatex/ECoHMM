-module(echo_supervisor).
-behaviour(supervisor).

%% API
-export([start_link/0]).

%% Supervisor callbacks
-export([init/1]).

-define(SERVER, ?MODULE).

%%====================================================================
%% API
%%====================================================================

start_link() ->
    supervisor:start_link({local, ?SERVER}, ?MODULE, []).

%%====================================================================
%% Supervisor callbacks
%%====================================================================

init([]) ->
    Dog =
        {dog,
         {echo_server, start_link, [dog, "Woof!", 3000]},
         permanent, 2000, worker, dynamic},

    HungryCat =
        {hungry_cat,
         {echo_server, start_link, [hungry_cat, "*meow*", 4000]},
         permanent, 2000, worker, dynamic},

    HappyCat =
        {happy_cat,
         {echo_server, start_link, [happy_cat, "*purrrrrr*", 1500]},
         permanent, 2000, worker, dynamic},

    Children = [Dog, HungryCat, HappyCat],
    {ok, {{one_for_one, 3, 10}, Children}}.

%%====================================================================
%% Internal functions
%%====================================================================